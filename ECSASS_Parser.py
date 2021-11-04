"""
ECSASS PARSER
(version 1.0)
by Angelo Chan

The Ebert-Chan Sequence Assembly Syntax System Parser module, implemented in
Python.

The ECSASS system is a compact system for storing genetic sequence data,
designed for scenarios where the same sequence is duplicated multiple tims,
and some copies may be altered differently to other copies. The ECSASS system
records the original sequences which the final products are based on, and any
changes or alterations which have been made to them, rather than a hard copy of
the altered sequences themselves.

The ECSASS system allows multiple "changes" to the sequences to be made with
ease, and be displayed in a human-readable format while taking up very little
hard disk space. The downside is that ECSASS sequences need to be assembled
first before being able to be used as a standard genetic sequence, which will
increase the runtime of software which bridges between the ECSASS format and the
direct basepair by basepair recroding format.



SYNTAX:
    
    FILE(<name>)
        Refer to the genetic sequence from a file with filename <name>.
    
    SEQ(<seq>)
        Directly provide a genetic sequence in plain text format.
    
    INV()
        Invert the sequence specified within the brackets.
    
    +
        Join or concatenate two sequences together.
    
    *<X>
        Duplicate a sequence <X> number of times.
    
    ~+
        Overlap-join or overlap-concatenate two sequences together. The last N
        nucleotides of the first sequence and the last N nucleotides of the last
        sequence which align will "merge" and only appear once. 
    
    ~*<X>
        Overlap-duplicate a sequence <X> number of times. The last N nucleotides
        of the sequence and the first X nucleotides which align will "merge" and
        only appear once. 

            N
                (In the case of ~+ and ~*<X>)
                N is either the largest or smallest window size within a
                specified range for which there is a valid overlap. Whether the
                largest or smallest window size is used is specified by the
                user, as is the minimum and maximum window sizes and the
                maximum permissible number of mismatches.

    [<X>:<Y>]
        Standard string slicing. Remove all nucleotides before index <X>, not
        inclusive, and after index <Y>, inclusive. If no number is provided
        before the colon, all nucleotides up to index <Y>, not inclusive, will
        be kept. If no number is provided after the colon, all nucleotides after
        index <X>, inclusive, will be kept.
        The first position is index 0. Negative numbers can be used to indicate
        the index when counting backwards from the end.

    ![<X>:<Y>]
        Truncation string slicing. All nucleotides after index <X>, inclusive,
        and before index <Y>, not inclusive. The first position is index 0.
        Negative numbers can be used to indicate the index when counting
        backwards from the end.

    ()
        Curved brackets can be used to enforce order of operations.



EXAMPLES:

    SEQ(NNNNNAACC)~+SEQ(AACCNNNNN)
        
        ->  NNNNNAACCNNNNN
    
    SEQ(NNNNNAACC)~+SEQ(AACCNNNNN)
        
        ->  NNNNNAACCNNNNN
                
    SEQ(AACCNNNNNAACC)~*3
        
        ->  AACCNNNNNAACCNNNNNAACCNNNNNAACC
                
    SEQ(AAACCCGGGTTT)[3:9]
        
        ->  CCCGGG
                
    SEQ(AAACCCGGGTTT)[3:]
        
        ->  CCCGGGTTT
                
    SEQ(AAACCCGGGTTT)[:-3]
        
        ->  AAACCCGGG
                
    SEQ(AAACCCGGGTTT)![-9:-3]
        
        ->  AAATTT
                
    SEQ(AAACCCGGGTTT)![-9:-3]
        
        ->  AAATTT
                
    INV(SEQ(AAACCCCCC))
        
        ->  GGGGGGTTT
"""

# Imported Modules #############################################################

import sys
import os

from NSeq_Match import *
from _Command_Line_Parser import *



# Sets #########################################################################

LIST__nuc = set(["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"])

LIST__fasta = set(["fa", "FA", "Fa", "fasta", "FASTA", "Fasta"])



# Functions ####################################################################

def Parse_ECSASS(ECSASS_seq, seq_folders, overhang_min, overhang_max, error_max,
            highest_preferred):
    """
    Parse and return an ECSASS sequence and transform it into a normal genetic
    sequence.
    
    @ECSASS_seq
            (str)
            The input ECSASS sequence to be parsed. Requires that none of the
            files referenced have whitespaces in their file names. If file names
            contain whitespaces, users may need to resort to _Parse_ECSASS(),
            which may be rendered nonfunctional by the presence of other
            whitespace characters.
    @input_sequences
            (list<str - dirpath>)
            The folders containing the FASTA files containing the sequences
            which will be used to generate the final sequences. In the event of
            a name conflict, the folders listed earlier in this list will be
            given priority.
    @overhang_min
            (int)
            The minimum window size allowed for overlap-joins and
            overlap-duplicates.
    @overhang_max
            (int)
            The maximum window size allowed for overlap-joins and
            overlap-duplicates.
    @error_max
            (int)
            The maximum number of mismatches permitted for two sequences to
            qualify for an overlap-join or an overlap-duplicate.
    @highest_preferred
            (bool)
            Whether or not the largest qualifying window size will be used for
            overlap-joins and overlap-duplicates. If False, the smallest
            qualifying window will be used instead.
    
    Parse_ECSASS(str, str, int, int, int, bool) -> str
    """
    ECSASS_seq = Strip_Whitespaces(ECSASS_seq)
    return _Parse_ECSASS(ECSASS_seq, seq_folders, overhang_min, overhang_max, error_max,
            highest_preferred)

def Strip_Whitespaces(string):
    """
    Return [string] with all whitespaces removed.
    """
    sb = ""
    for c in string:
        if c != " ": sb += c
    return sb

def _Parse_ECSASS(ECSASS_seq, seq_folders, overhang_min, overhang_max, error_max,
            highest_preferred):
    """
    The recursive, subfunction of Parse_ECSASS(). Assumes [ECSASS_seq] contains
    no whitespace characters.
    """
    # Flags and counters
    parts = []
    operations = []
    sb = ""
    bracket_level = 0
    # Parsing infrastructure
    length = len(ECSASS_seq)
    # Traverse string
    i = 0
    while i < length:
        char = ECSASS_seq[i]
        # Brackets
        if char == "(":
            bracket_level += 1
        elif char == ")":
            bracket_level -= 1
            if bracket_level < 0:
                raise Exception("ERROR: Unmatched closing bracket.")
            elif bracket_level == 0:
                sb = _Parse_ECSASS(sb[1:], seq_folders, overhang_min,
                        overhang_max, error_max, highest_preferred)
                parts.append(sb)
                sb = ""
        # In brackets
        if bracket_level:
            sb += char
        # File
        elif ECSASS_seq[i:i+5] == "FILE(":
            i += 5
            char = ECSASS_seq[i]
            while char and char != ")":
                sb += char
                i += 1
                char = ECSASS_seq[i]
            seq = Get_Seq_From_File(sb, seq_folders)
            parts.append(seq)
            sb = ""
        # Sequence
        elif ECSASS_seq[i:i+4] == "SEQ(":
            i += 4
            char = ECSASS_seq[i]
            while char and char != ")":
                # if char in LIST__nuc: raise Exception
                sb += char
                i += 1
                char = ECSASS_seq[i]
            parts.append(sb)
            sb = ""
        # Invert
        elif ECSASS_seq[i:i+4] == "INV(":
            i += 4
            bracket_level += 1
            while bracket_level > 0:
                try:
                    char = ECSASS_seq[i]
                except:
                    raise Exception("ERROR: Unmatched closing bracket.")
                # Brackets
                if char == "(":
                    bracket_level += 1
                elif char == ")":
                    bracket_level -= 1
                # Other
                sb += char
                i += 1
            sb = sb[:-1] # Remove trailing bracket
            sb = _Parse_ECSASS(sb, seq_folders, overhang_min, overhang_max,
                    error_max, highest_preferred)
            sb = Get_Complement(sb, True)
            parts.append(sb)
            sb = ""
        # Prep for next loop iteration
        i += 1
    # Unclosed loop
    if bracket_level:
        raise Exception("ERROR: Unmatched opening bracket.")
    # Combine Parts
    while operations:
        operation = operations[0]
        pass
    # Return
    return parts[0]

def Get_Seq_From_File(name, folders):
    """
    Get the genetic sequence from a file with [name] from folders. The folders
    will be scanned in order of listing, with the folders listed earlier
    effectively being given priority in the event of a naming clash.
    
    @name
            (str)
            The name of the file, such as a gene name.
    @folders
            (list<str - dirpath>)
            A list of folder paths, which will be scanned.
    
    Get_Seq_From_File(str, list<str>) -> str
    """
    # Attempt to locate file
    final_path = ""
    for folder in folders:
        file_names = os.listdir(folder)
        for file_name in file_names:
            file_name_, ext = Get_File_Name_Ext(file_name)
            if file_name_ == name and ext in LIST__fasta:
                searching = False
                if not final_path: final_path = folder + "//" + file_name
    # File not located
    if not final_path: raise Exception("No file named \"{S}\" found.".format(
            S = name))
    # File located, extract sequence
    seq = ""
    f = open(final_path, "U")
    f.readline()
    for line in f:
        if line and (line[-1] == "\n" or line[-1] == "\r"): line = line[:-1]
        seq += line
    f.close()
    return seq


