#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#
# simple field mapper, provide mapping files, supports chaining of mappings
#

__author__ = "johdro, weimann"

from sys import stderr, stdout, exit, argv, stdin
import re

#TODO add non-interactive mode 

class SequentialMapper:


    def __init__(self, fh, mfield=0, is_incomplete_mapping=False, has_duplicates=False):
        self._fh = fh
        self._is_incomplete_mapping = is_incomplete_mapping
        self._mfield = mfield
        self._lastkey = None
        self._curkey, self._curval = self.next()
        self._matched = False
        self._has_duplicates = has_duplicates

    def __getitem__(self, key):
        if key is  None:
            self._curkey, self_curval = self.next()
            return self._curval
        if not self.an(key, self._curkey) and self.an(self._curkey, key) and self._lastkey is not None and self._lastkey == key:
            if self._has_duplicates:
                return self._curval
            else:
                stderr.write("%s is duplicated and allow duplicate option (-d is not put)"%key)
                raise ValueError
        while self.an(key, self._curkey):
            try:
                # check if cur_key has already been matched, otherwise output
                # curkey as it couldn't be mapped
                if not self._matched:
                    print "NA\t%s" % self._curkey
                self._matched = False
                self._curkey, self._curval = self.next()
                # check if after update of current key, key still can be
                # matched, output NA if not
                if self.an( self._curkey, key) and self._is_incomplete_mapping:
                    return ["NA"]
                elif not self._is_incomplete_mapping:
                    stderr.write("%s could not be mapped" % key)
                    raise ValueError
            except StopIteration:
                if not self._is_incomplete_mapping:
                    raise ValueError
                else: return ["NA"]
        # check if the key still can be matched, output NA if not
        if self.an(self._curkey, key) and self._is_incomplete_mapping:
            return ["NA"]
        elif not self._is_incomplete_mapping:
            stderr.write("%s could not be mapped" % key)
            raise ValueError
        # key and current key matched
        self._matched = True
        self._lastkey = key
        return self._curval

    @staticmethod
    def an(key1, key2):
        """replace all non alphanumerical characters with underlines and make keys the same size to compare"""
        #print key1, key2
        if len(key1) > len(key2): key2+="_"*(len(key1) - len(key2))
        else: key1+="_"*(len(key2) - len(key1))
        #print re.sub("[^A-Za-z0-9]", "", key1).lower() > re.sub("[^A-Za-z0-9]", "", key2).lower()
        return re.sub("[^A-Za-z0-9]", "", key1).lower() > re.sub("[^A-Za-z0-9]", "", key2).lower()

    def next(self):
        fh_list = self._fh.next().rstrip("\n").split("\t")
        # return tuple of key (first element) and value (a list of all other
        # elements)
        return tuple([fh_list[self._mfield], fh_list[0: len(fh_list)]])


def transMapping(value, mappers):  # like foldl()
    for m in mappers:
        value = m[value]
    return value


def usage():
    stderr.write(
        "Usage: %s -f --field (default 1) -m mapping_field (default 1) -s --mapping-sorted 1.tsv  -d --has_duplicates (default false) -i --incomplete_mapping (default false) < 2.tsv\n" %
        argv[0])
    stderr.write(
            """The input and the mapping have to be sorted by GNU sort -d (dictionary sort) for this script to work, the mapping has to be additionally sorted uniquly by sort -u, if a mapping key is not matched by an input key it will be outputed as \"NA\tmapping_key\"\n All non-key fields of the mapping and the input file will be appended to the output. Use cut or awk to further restrict your output.\n Note: All non-alphanumerical characters in the key column are ignored by GNU sort and the python compare method!!\n""")
    stderr.write("If the incomplete mapping option is given and a mapping key cannot be matched this script will output it as  \"input_key\tNA\"\n")


if __name__ == "__main__":

    import getopt
    #s=an("test#*123_")
    #assert(s=="test__123_")
    # handle broken pipes
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    # parse command line options
    try:
        opts, args = getopt.getopt(
            argv[
                1:], "hf:s:idm:", [
                "help", "field", "mapping-sorted=", "incomplete-mapping", "has_duplicates","mapping_field"])
        #opts, args = getopt.getopt(argv[1:], "hf:s:m:i", ["help", "field", "mapping-sorted=", "mapping="])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        exit(2)
    if len(argv) == 1:
        usage()
        exit(2)
    # default parameters
    mappers = []
    ifield_num =0
    mfield_num = 0
    is_incomplete_mapping = False
    has_duplicates = False
    for o, a in opts:
        if o in ("-i", "--incomplete-mapping"):
            is_incomplete_mapping = True
        elif o in ("-d", "--has_duplicates"):
            has_duplicates = True
        elif o in ("-m", "--mapping_field"):
            mfield_num = int(a) - 1

    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            exit()
        elif o in ("-f", "--input_field"):
            ifield_num = int(a) - 1
        elif o in ("-s", "--mapping-sorted"):
            mappers.append(
                SequentialMapper(
                    open(a,"r"), has_duplicates = has_duplicates, mfield=mfield_num, is_incomplete_mapping = is_incomplete_mapping
                )
            )
        # elif o in ("-m", "--mapping"):
        #    mappers.append(dict([line.rstrip("\n").split("\t", 2)[:2] for line in open(a, "r")]))
        elif o in ("-m", "--field"):
            pass
        elif o in ("-i", "--incomplete-mapping"):
            pass
        elif o in ("-d", "--has_duplicates"):
            pass
        else:
            assert False, "unhandled option"

    for line in stdin:
        if line and line[0] == "#":
            continue
        fields = line.rstrip('\n').split("\t")
        fields= fields[0:10]
        fields = fields + transMapping(fields[ifield_num], mappers)
        stdout.write("\t".join(fields))
        stdout.write("\n")
    if is_incomplete_mapping:
        while True:
            try:
                stdout.write("NA\t%s\n"%"\t".join(mappers[0][None]))
            except StopIteration:
                exit(0)
