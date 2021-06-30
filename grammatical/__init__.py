from re import finditer, DOTALL, MULTILINE
import json
from functools import partial
import logging


def mapping_depth_first(self, codons, order=None):
    """

    :param self:
    :param codons:
    :param order:
    :return:
    """
    pro = [self.start]
    nt = [0]

    if order is None:
        order = range(len(codons))

    i = 0
    while i < len(order) and len(nt) > 0:
        codon = int(codons[order[i]])
        pos = codon % len(self.productions[pro[nt[0]]])
        ipos = nt[0] + 1
        for item in self.productions[pro[nt[0]]][pos]:
            pro.insert(ipos, item)
            ipos += 1

        del pro[nt[0]]
        nt = [x for (x, item) in enumerate(pro) if item in self.non_terminals]
        i += 1

    if len(nt) != 0:
        return None

    #return "".join(pro)
    
    npro = []
    count = self.start_count
    for item in pro:
        if self.replace_chart in item:
            npro.append(item.replace(self.replace_chart, str(count)))
            count += 1
        else:
            npro.append(item)
    return "".join(npro)


def mapping_depth_first2(self, codons, order=None):
    """

    :param self:
    :param codons:
    :param order:
    :return:
    """

    pro = [self.start]
    position = 0
    if order is None:
        order = range(len(codons))

    for order_codon in order:
        codon = codons[order_codon]
        for i in range(position, len(pro)):
            if pro[i] in self.non_terminals:
                pos = codon % len(self.productions[pro[i]])
                ipos = i + 1
                for item in self.productions[pro[i]][pos]:
                    pro.insert(ipos, item)
                    ipos += 1

                del pro[i]
                position = i
                break

    for item in pro:
        if item in self.non_terminals:
            return None

    return "".join(pro)


def mapping_breadth_first(self, codons, order=None):
    """

    :param self:
    :param codons:
    :param order:
    :return:
    """

    pro = [self.start]
    i = 0
    if order is None:
        order = range(len(codons))

    for order_codon in order:
        codon = int(codons[order_codon])

        step = 0
        while True:
            if i >= len(pro):
                i = 0
                step += 1
            if step > 1:
                break

            if pro[i] in self.non_terminals:
                pos = codon % len(self.productions[pro[i]])
                ipos = i
                for item in self.productions[pro[i]][pos]:
                    ipos += 1
                    pro.insert(ipos, item)

                del pro[i]
                i = ipos
                break
            i += 1

        if step > 1:
            break

    for item in pro:
        if item in self.non_terminals:
            return None

    return "".join(pro)


def mapping_pigrammatical(self, codons, order=None):
    """

    :param self:
    :param codons:
    :param order:
    :return:
    """

    pro = [self.start]
    nt = [0]

    size = int(len(codons) / 2)
    if order is None:
        order = range(size)

    list_nt, list_codons = codons[0::2], codons[1::2]

    i = 0
    while i < size and len(nt) > 0:
        nt_value = int(list_nt[order[i]])
        ntpos = nt_value % len(nt)
        codon = int(list_codons[order[i]])
        pos = codon % len(self.productions[pro[nt[ntpos]]])
        ipos = nt[ntpos] + 1

        for item in self.productions[pro[nt[ntpos]]][pos]:
            pro.insert(ipos, item)
            ipos += 1

        del pro[nt[ntpos]]

        nt = [x for (x, item) in enumerate(pro) if item in self.non_terminals]
        i += 1

    if len(nt) != 0:
        return None

    return "".join(pro)


class Grammar:
    """
    """
    productions = {}
    terminals = []
    non_terminals = []
    start = None

    def __init__(self, mapping=mapping_depth_first, file=None, grammar=None, start_count=0, replace_chart="n"):
        # setattr(self.__class__, "mapping", mapping)
        self.mapping = partial(mapping, self)
        self.start_count = start_count
        self.replace_chart = replace_chart
        if grammar is not None:
            self.productions = grammar.productions
            self.terminals = grammar.terminals
            self.non_terminals = grammar.non_terminals
            self.start = grammar.start

        elif file is not None:
            ruleregex = '(?P<rulename><\S+>)\s*::=\s*(?P<production>(?:(?=\#)\#[^\r\n]*|(?!<\S+>\s*::=).+?)+)'
            productionregex = '(?=\#)(?:\#.*$)|(?!\#)\s*(?P<production>(?:[^\'\"\|\#]+|\'.*?\'|".*?")+)'
            productionpartsregex = '\ *([\r\n]+)\ *|([^\'"<\r\n]+)|\'(.*?)\'|"(.*?)"|(?P<subrule><[^>|\s]+>)|([<]+)'

            if file.endswith(".json"):
                with open(file, 'r') as fp:
                    gramm = json.load(fp)
                    self.start = gramm["start"]
                    self.terminals = gramm["terminals"]
                    self.non_terminals = gramm["non_terminals"]
                    for key in gramm["productions"]:
                        if key in self.productions:
                            continue

                        pro = []
                        for p in gramm["productions"][key]:
                            subpro = []
                            for sub_p in finditer(productionpartsregex, p.strip()):
                                if sub_p.group('subrule') is not None:
                                    subpro.append(sub_p.group('subrule'))
                                else:
                                    sub = ''.join(
                                        [part.encode().decode('unicode-escape') for part in sub_p.groups() if
                                         part]).strip()
                                    if sub != '':
                                        subpro.append(sub)
                            pro.append(subpro)
                        self.productions[key] = pro

            elif file.endswith(".bnf"):
                with open(file, 'r') as fp:
                    for line in fp:
                        for rule in finditer(ruleregex, line, DOTALL):
                            if self.start is None:
                                self.start = rule.group('rulename')

                            if rule.group('rulename') not in self.non_terminals:
                                self.non_terminals.append(rule.group('rulename'))

                            if rule.group('rulename') in self.productions:
                                continue

                            pro = []
                            for p in finditer(productionregex, rule.group('production'), MULTILINE):
                                if p.group('production') is None:
                                    continue

                                subpro = []
                                for sub_p in finditer(productionpartsregex, p.group('production').strip()):
                                    if sub_p.group('subrule') is not None:
                                        subpro.append(sub_p.group('subrule'))
                                    else:
                                        sub = ''.join(
                                            [part.encode().decode('unicode-escape') for part in sub_p.groups() if
                                             part]).strip()
                                        if sub != '':
                                            subpro.append(sub)
                                pro.append(subpro)
                            self.productions[rule.group('rulename')] = pro
            else:
                raise Exception("Is neccesary to give a file or a Grammar instance")

    def __str__(self):
        return "Start ::= {}\nTerminals ::= {}\nNon Terminals ::= {}\nProduction ::={}".format(self.start, json.dumps(self.terminals, sort_keys=True, indent=4 * ' '), json.dumps(self.non_terminals, sort_keys=True, indent=4 * ' '), json.dumps(self.productions, sort_keys=True, indent=4 * ' '))

    def mapping(self, codons, order=None):
        pass
