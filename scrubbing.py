# -*- coding: utf-8 -*-
# @Author: Martin Grunnill
# @Date:   2023-12-15 12:56:21
# @Last Modified by:   Martin Grunnill
# @Last Modified time: 2023-12-15 15:28:51
"""
Date Created: April 25th 2023
Author: Martin Grunnill 
"""
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import OrderedDict, defaultdict
import re
from datetime_funcs import is_date, date_to_year_fraction, string_convert_date_format
from dateutil.parser import parse
import pandas as pd


def write_seq_dict_to_file(sequence_dict, file_name):
    """_summary_

    Parameters
    ----------
    sequence_dict : dict
        Dictionary of genome/protiome sequences. {tag/ID: sequence}
    file_name : str
        Path to write files to.
    """
    with open(file_name, "w+") as output_file:
        for id, sequence in sequence_dict.items():
            output_file.write('>'+id + "\n" + sequence + "\n")

def fasta_to_dict(input_file, delimiter='|'):
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if delimiter is not None:
            tag = str(seq_record.description).split(delimiter,1)[0]
        else:
            tag = str(seq_record.description)
        sequence_dict[tag] = str(seq_record.seq)
    
    return sequence_dict


def fasta_to_Bankit_format(input_file, output_file, metadata):
    """Convert fasta file to Bankit format.

    See following for Bankit format guidelines:
        https://www.ncbi.nlm.nih.gov/WebSub/html/help/fasta.html and
        https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html#modifiers


    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    metadata : pandas.DataFrame
        Sequence tag must correspond to index of dataframe. To set index use metadata.set_index('ID',inplace=True).
    """
    seq_dict = fasta_to_dict(input_file)
    conv_seq_dict = {}
    for count, (tag, seq) in enumerate(seq_dict.items(),1):
        metadata_entry = metadata.loc[tag].to_dict()
        metadata_entry = ['['+column+'='+row+']' for column, row in metadata_entry.items()]
        metadata_entry = ' '.join(metadata_entry)
        conv_seq_dict['Seq'+str(count)+' '+metadata_entry] = seq
    
    write_seq_dict_to_file(conv_seq_dict, output_file)

def seq_clean(input_file, output_file, accept = 'iupac_codes', replacement='n'):
    """Convert unacceptable characters from sequences in fasta file entries to a replacement character.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    accept : str, optional
        Nucleotides or iupac_codes (plus '-'), by default 'iupac_codes'.
    replacement : str, optional
        Replacement character, by default 'n'.

    Raises
    ------
    AssertionError
        _description_
    """
    nucleotides = 'g,c,a,t,u'
    iupac_codes = nucleotides + 'r,y,m,k,s,w,h,b,v,d,n,-'
    iupac_without_alignment =  nucleotides + 'r,y,m,k,s,w,h,b,v,d,n'
    if accept == 'nucleotides':
        accepted = nucleotides
    elif accept == 'iupac_codes':
        accepted = iupac_codes
    elif accept == 'iupac_without_alignment':
        accept = iupac_without_alignment
    else:
        raise AssertionError('Argument accept only takes the values of "nucleotides", "iupac_codes" or "iupac_without_alignment"')
    accepted = '[^.' + accepted + ']'
    
    sequence_dict = OrderedDict()

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        # Clean sequence and place in dictionary
        sequence_dict[str(seq_record.description)] = re.sub(accepted, replacement, str(seq_record.seq).lower())
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def seq_filter(input_file, output_file, min_length=0, por_n=100):
    """Filter fasta file entries by quality.

    Adapted from https://biopython.org/wiki/Sequence_Cleaner

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    min_length : int, optional
        Minimum length, by default 0
    por_n : number, optional
        _Percentage of characters in n acceptable, by default 100
    """
    # Create our hash table to add the sequences
    sequence_dict = OrderedDict()

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq)
        # Check if the current sequence is according to the user parameters
        if (
            len(sequence) >= min_length
            and (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n
        ):
            sequence_dict[seq_record.description] = sequence

    # Write the filtered sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def tag_convert_delimiter(input_file, output_file, delimiter='_', replacement='|'):
    """Convert delimiter in fasta file entry tags (IDs) to replacement character.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    delimiter : str, optional
        Delimiter between rest of tag (ID) and date, by default '_'.
    replacement : str, optional
        Replacement delimiter by default '/'
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).replace(delimiter, replacement)
        # Place in dictionary
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def tag_convert_date_delimiter(input_file, output_file, last_delimiter='_', replacement='|'):
    """Convert date delimiter in fasta file entry tags (IDs) to replacement character.

    Assumes date is at end of tags (IDs).

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    last_delimiter : str, optional
        Delimiter between rest of tag (ID) and date, by default '_'.
    replacement : str, optional
        Replacement delimiter by default '/'
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).rsplit(last_delimiter,1)
        if is_date(tag[-1]):
            tag = replacement.join(tag)
        else:
            tag = str(seq_record.description)
        # Place in dictionary
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)

def tag_convert_last_delimiter(input_file, output_file, last_delimiter='_', replacement='|'):
    """Convert last delimiter in fasta file entry tags (IDs) to replacement character.


    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    last_delimiter : str, optional
        Delimiter between rest of tag (ID) and date, by default '_'.
    replacement : str, optional
        Replacement delimiter by default '/'
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).rsplit(last_delimiter,1)
        tag = replacement.join(tag)
        # Place in dictionary
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def remove_if_no_date(input_file, output_file, last_delimiter='|'):
    """Remove fata file entry if no date present at end tag (ID).

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    last_delimiter : str, optional
        Delimiter, by default '|'
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).rsplit(last_delimiter,1)
        if is_date(tag[-1]):
            sequence_dict[str(seq_record.description)] = str(seq_record.seq)

    # Write the filtered sequences
    write_seq_dict_to_file(sequence_dict, output_file)

def tag_convert_date_format(input_file, output_file, last_delimiter='|', old_format='%m-%d-%Y', new_format='%Y-%m-%d'):
    """Convert date format of fasta entery tag, assuming date is at end of tag.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    last_delimiter : str, optional
        Delimiter, by default '|'
    old_format : str, optional
        Old date format, by default '%m-%d-%Y'
    new_format : str, optional
        New date format, by default '%Y-%m-%d'
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).rsplit(last_delimiter,1)
        date = tag[-1]
        tag = tag[0]
        date = string_convert_date_format(date, old_format=old_format, new_format=new_format)
        tag += last_delimiter + date     
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)



def tag_date_to_year_fraction(input_file, output_file, last_delimiter='|', yearfirst=True):
    """Change fasta file entry tags (IDs) date section to year fraction (decimal) format.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    last_delimiter : str, optional
        Character used for last_delmiter, by default '/'.

    Raises
    ------
    ValueError
        Raised if tag has 4 charachters after delimiter and they are not digits.
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).rsplit(last_delimiter,1)
        date = tag[-1]
        tag = tag[0]
        if len(date)==4:
            if not date.isdigit():
                raise ValueError('If ' + tag + last_delimiter + date + ' least segment (date) is 4 characters they must be digits')
            tag += last_delimiter + date + '.5'
        else:
            try:
                date = parse(date, yearfirst=yearfirst)
            except:
                raise ValueError(str(seq_record.id) + ' date ('+date+') is not recognised as a date format.')
            year_fraction = date_to_year_fraction(date)
            tag += last_delimiter + str(round(year_fraction, 3))       

        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)

def tag_strip_after_delimiter(input_file, output_file, delimiter='|'):
    """Strip fasta file entry tags (IDs) after delimiter.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    delimiter : str, optional
        Character used for delmiter, by default '|'.
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).split(delimiter,1)[0]
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def remove_if_accession_in_list(input_file, output_file, accession_list, delimiter='|'):
    """Remover if accession number is in list, assumes accession number is first in fasta entry tag.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    accession_list : list of strings or string
        List of accession numbers or single accession number.
    delimiter : str, optional
        Character used for delmiter, by default '|'.
    """
    if isinstance(accession_list,str):
        accession_list = [accession_list]
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.id)
        first_part_tag = tag.split(delimiter,1)[0]
        if first_part_tag not in accession_list: 
            sequence_dict[str(seq_record.description)] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)

def remove_if_accession_not_in_list(input_file, output_file, accession_list, delimiter='|'):
    """Remover if accession number is not in list, assumes accession number is first in fasta entry tag.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    accession_list : list of strings or string
        List of accession numbers or single accession number.
    delimiter : str, optional
        Character used for delmiter, by default '|'.
    """
    if isinstance(accession_list,str):
        accession_list = [accession_list]
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.id)
        first_part_tag = tag.split(delimiter,1)[0]
        if first_part_tag in accession_list: 
            sequence_dict[str(seq_record.description)] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def keep_if_tag_startswith(input_file, output_file, string):
    """Keep fasta entry if the associated tag startwith string.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    string : str
        String tags should start with.
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.id)
        if tag.startswith(string): 
            sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)


def remove_if_tag_startswith(input_file, output_file, string):
    """Remove fasta entry if the associated tag startwith string.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    string : str
        String tags should start with.
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.id)
        if not tag.startswith(string): 
            sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)

    
def move_entry_to_start(input_file, output_file, string):
    """Move fasta entry to beginning of fasta file if the tag starts with string.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    string : str
        String tags should start with.

    Raises
    ------
    ValueError
        _description_
    """
    description_searching_for = None
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.id)
        if tag.startswith(string):
            description_searching_for = str(seq_record.description)
        sequence_dict[str(seq_record.description)] = str(seq_record.seq)
    
    if description_searching_for is not None:
        sequence_dict.move_to_end(description_searching_for, last=False)
    else:
        raise ValueError('Entry for ' + accession_number + ' not found in ' +input_file)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)
            
    
def add_date_string_to_tag(input_file, output_file, dates_dict, delimiter='|'):
    """Add date(s) to end of fasta entry tag(s). 

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    dates_dict : dict
        Dictionary {accession_number: date}
    delimiter : str, optional
        Character used for delmiter, by default '|'.
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        accession_number = str(seq_record.id).split(delimiter,1)[0]
        tag = str(seq_record.description)
        if accession_number in dates_dict:
            tag += dates_dict[accession_number]
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)



def add_info_string_to_tag_pen_place(input_file, output_file, info, delimiter='|'):
    """Add information in a new penultimate field of fasta entry tags.

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    info : str or dict
        If string is given information is added to all fasta enties. 
        If dictionary the dictionary is used to map the information via accesion numbers.
    delimiter : str, optional
        Character used for delmiter, by default '|'.
    """
    sequence_dict = OrderedDict()
    if isinstance(info, dict):
        info_is_dict = True
    else:
        info_is_dict = False
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        accession_number = str(seq_record.id).split(delimiter,1)[0]
        tag = str(seq_record.description)
        if info_is_dict:
            if accession_number in info:
                tag_split = tag.rsplit(delimiter,1)
                tag = tag_split[0] + delimiter + info[accession_number] + delimiter + tag_split[1]
        else:
            tag_split = tag.rsplit(delimiter,1)
            tag = tag_split[0] + delimiter + info + delimiter + tag_split[1]
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the ccomverted sequences
    write_seq_dict_to_file(sequence_dict, output_file)

def sequence_len_df(input_file):
    """Report sequence lengths of fasta entires.

    Parameters
    ----------
    input_file : str
        Path to input file.

    Returns
    -------
    pandas.DataFrame
        Dataframe of accession numbers and sequence lengths.
    """
    df_contents_dict = {'identifiers': [],
                    'lengths' : []}
    with open(input_file) as fasta_file:  # Will close handle cleanly
        for title, sequence in SimpleFastaParser(fasta_file):
            df_contents_dict['identifiers'].append(title.split(None, 1)[0])  # First word is ID
            df_contents_dict['lengths'].append(len(sequence))
    
    
    return pd.DataFrame.from_dict(df_contents_dict)

def search_for_clad_strings(input_file, clade_srings, delimiter='_'):
    """_summary_

    Parameters
    ----------
    input_file : _type_
        _description_
    clade_srings : _type_
        _description_
    delimiter : str, optional
        _description_, by default '_'

    Returns
    -------
    _type_
        _description_
    """
    clade_dict = {clade_string: [] for clade_string in clade_srings}
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.id)
        first_part_tag = tag.split(delimiter,1)[0]
        for clade_string in clade_srings:
            clade_string_with_delim = delimiter+clade_string+delimiter
            if clade_string_with_delim in str(seq_record.description):
                clade_dict[clade_string].append(first_part_tag)
    
    return clade_dict


def replace_with_sequence_from(input_file, replacements_file, output_file,
                               input_delimiter, replacement_delimiter, replacement_tag_suffix='',
                               replace_if_unknown_symbol=False):
    nucleotides = 'g,c,a,t,u'
    iupac_codes = nucleotides + 'r,y,m,k,s,w,h,b,v,d,n,-'
    replacement_dict = {}
    for seq_record in SeqIO.parse(replacements_file, "fasta"):
        tag = str(seq_record.id)
        first_part_tag = tag.split(replacement_delimiter,1)[0]
        sequence = str(seq_record.seq)
        if not replace_if_unknown_symbol:
            if all((charecter in iupac_codes) for charecter in sequence):
                replacement_dict[first_part_tag+replacement_tag_suffix] = sequence
        else:    
            replacement_dict[first_part_tag+replacement_tag_suffix] = sequence
    
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        accession_number = str(seq_record.id).split(input_delimiter,1)[0]
        tag = str(seq_record.description)
        if accession_number in replacement_dict:
            sequence_dict[tag] = replacement_dict[accession_number]
        else:
            sequence_dict[tag] = str(seq_record.seq)
    
    write_seq_dict_to_file(sequence_dict, output_file)


def tag_add_suffix_to_first_part(input_file, output_file, delimiter='|', suffix='.1 '):
    """_summary_

    Parameters
    ----------
    input_file : str
        Path to input file.
    output_file : str
        Path to wrtie output file to.
    delimiter : str, optional
        Delimiter between rest of tag (ID) and date, by default '|'.
    suffix : str, optional
        Suffix to add to first part of tag, by default '.1'
    """
    sequence_dict = OrderedDict()
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).split(delimiter, 1)
        tag = tag[0] + suffix + delimiter + tag[1]
        # Place in dictionary
        sequence_dict[tag] = str(seq_record.seq)
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)   


def add_fasta_entries_if_not_in(input_file, to_be_added_file, output_file, delimiter='|'):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        accession = str(seq_record.description).split(delimiter, 1)[0]
        # Place in dictionary
        sequence_dict[accession] = (str(seq_record.description), str(seq_record.seq))
    

    for seq_record in SeqIO.parse(to_be_added_file, "fasta"):
        accession = str(seq_record.description).split(delimiter, 1)[0]
        if accession not in sequence_dict:
            sequence_dict[accession] = (str(seq_record.description), str(seq_record.seq))
    
    sequence_dict = OrderedDict([entry for entry in sequence_dict.values()])
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file) 


def ids_per_sequence(input_file, delimiter='|'):
    """Get ids per sequence that is the same (identical) in fasta file.

    Parameters
    ----------
    input_file : string
        Path to fasta file.
    delimiter : str, optional
        Delimiter between rest of tag (ID) and date, by default '|'.

    Returns
    -------
    defaultdict of lists
        Keys are sequence values are a list of ids.

    """    
    sequence_dict = defaultdict(list)
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag =str(seq_record.description).split(delimiter, 1)[0]
        sequence_dict[str(seq_record)].append(tag)
    return sequence_dict
    

def get_ids_of_duplicated_sequences(input_file, delimiter='|', return_full_dict=False):
    """Retrieves ids of sequences that are duplicated in a fasta file.

    Parameters
    ----------
    input_file : string
        Path to fasta file.
    delimiter : str, optional
        _description_, by default '|'
    return_full_dict : bool, optional
        Modifies return from dictionary {sequences: list of ids} or a list of sub lists of ids, by default False

    Returns
    -------
    if return_full_dict:
        dict of lists
            Keys are sequence values are a list of ids.
    else:
        list of lists:
            Sublist enties are ids.
    """
    sequence_dict = ids_per_sequence(input_file, delimiter)
    if return_full_dict:
        return {seq:id_list for seq, id_list in sequence_dict.items() if len(id_list)>1}
    else:
        return [id_list for id_list in sequence_dict.values() if len(id_list)>1]

        

def remove_if_any_substrings_in_tag(input_file, output_file, sustrings, delimiter='|',
                                    provid_accessions_of_removed=True):
    sequence_dict = OrderedDict()
    if provid_accessions_of_removed:
        removed = []
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        if not any(substring in tag for substring in sustrings):
            sequence_dict[tag] = str(seq_record.seq)
        elif provid_accessions_of_removed:
            accession = str(seq_record.description).split(delimiter, 1)[0]
            removed.append(accession)
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)
    if provid_accessions_of_removed:
        return removed


def accession_comparison_two_fastas(fasta1,fasta2,delimiter='|'):
    accessions_in_fasta1 = {str(seq_record.description).split(delimiter, 1)[0]
                            for seq_record in SeqIO.parse(fasta1, "fasta")}
    accessions_in_fasta2 = {str(seq_record.description).split(delimiter, 1)[0]
                            for seq_record in SeqIO.parse(fasta2, "fasta")}
    comparison_dict = {'intersection': accessions_in_fasta1.intersection(accessions_in_fasta2),
                       'Only in '+ fasta1: accessions_in_fasta1.difference(accessions_in_fasta2),
                       'Only in '+ fasta2: accessions_in_fasta2.difference(accessions_in_fasta1)}
    return comparison_dict

def tag_replace_substrings(input_file, output_file, replacement_dict):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        for old, new in replacement_dict.items():
            tag = tag.replace(old, new)
        sequence_dict[tag] = str(seq_record.seq)
    
    write_seq_dict_to_file(sequence_dict, output_file)

def seq_replace_substrings(input_file, output_file, replacement_dict):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        sequence = str(seq_record.seq)
        for old, new in replacement_dict.items():
            sequence = sequence.replace(old, new)
        sequence_dict[str(seq_record.description)] = sequence
    
    write_seq_dict_to_file(sequence_dict, output_file)

def remove_if_character_prop_over(input_file, output_file, prop, character, delimiter='|',
                                    provid_accessions_of_removed=True):
    sequence_dict = OrderedDict()
    if provid_accessions_of_removed:
        removed = []
    for seq_record in SeqIO.parse(input_file, "fasta"):
        sequence = str(seq_record.seq)
        char_prop = sequence.count(character)/len(sequence)
        if char_prop <= prop:
            sequence_dict[str(seq_record.description)] = sequence
        elif provid_accessions_of_removed:
            accession = str(seq_record.description).split(delimiter, 1)[0]
            removed.append(accession)
    
    # Write the cleaned sequences
    write_seq_dict_to_file(sequence_dict, output_file)
    if provid_accessions_of_removed:
        return removed


def tag_metadata_to_pandas(input_file, fields, delimiter='|'):
    records = []
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description).split(delimiter)
        records.append({fields[index]: item for index, item in enumerate(tag)})
    

    return pd.DataFrame.from_records(records)


def select_part_of_sequence(input_file, output_file, selection_dict, delimiter='|'):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        accession = tag.split(delimiter, 1)[0]
        if accession in selection_dict:
            start, end = selection_dict[accession]
            sequence = str(seq_record.seq)
            sequence_dict[tag] = sequence[start-1:end]
    
    write_seq_dict_to_file(sequence_dict, output_file)

def replace_tag(input_file, output_file, replacement_dict, delimiter='|'):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        accession = tag.split(delimiter, 1)[0]
        if accession in replacement_dict:
            sequence_dict[replacement_dict[accession]] = str(seq_record.seq)
        else:
            sequence_dict[tag] = str(seq_record.seq)
            
    write_seq_dict_to_file(sequence_dict, output_file)

def trim_aligned_to_match_length_of_sequence(input_file, output_file, selected_accession, delimiter='|'):
    found = False
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        accession = tag.split(delimiter, 1)[0]
        if accession == selected_accession:
            if found == True:
                raise ValueError('selected_accession appears more than once in fasta file.')
            found = True
            selected_seq = str(seq_record.seq)
            nucleotides = selected_seq.strip('-')
            gaps = selected_seq.split(nucleotides)
            if len(gaps)!=2:
                raise AssertionError('Should only be two gap regions.')
            start_index = len(gaps[0])
            end_index = start_index + len(nucleotides)
    
    if not found:
        raise ValueError('Accession ('+ selected_accession + ') not found')
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        seq = str(seq_record.seq)
        sequence_dict[tag] = seq[start_index:end_index]
            
    write_seq_dict_to_file(sequence_dict, output_file)

def translate_sequnces(input_file, output_file):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        seq = str(seq_record.seq.translate())
            
    write_seq_dict_to_file(sequence_dict, output_file)

def trim_sequences(input_file, output_file, start_index, end_index):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        seq = str(seq_record.seq)
        sequence_dict[tag] = seq[start_index:end_index]
            
    write_seq_dict_to_file(sequence_dict, output_file)


def tag_replace_substrings(input_file, output_file, replacement_dict):
    sequence_dict = OrderedDict()
    for seq_record in SeqIO.parse(input_file, "fasta"):
        tag = str(seq_record.description)
        for old, new in replacement_dict.items():
            tag = tag.replace(old, new)
        sequence_dict[tag] = str(seq_record.seq)
    
    write_seq_dict_to_file(sequence_dict, output_file)

# def merge_meta_data_fields_on_fasta(input_fastas, meta_data,
#                                     output_file,
#                                     fields=['country_original','collection_date'],
#                                     delimiter='/'):
#     sequence_dict = OrderedDict()
#     for input_file in input_fastas:
#         # Using the Biopython fasta parse we can read our fasta input
#         for seq_record in SeqIO.parse(input_file, "fasta"):
#             tag = str(seq_record.id)
#             field_data_list = [meta_data.loc(tag,field) for field in fields]
#             new_tag = [tag] + field_data_list
#             new_tag = delimiter.join(new_tag)
#             sequence_dict[new_tag] = str(seq_record.seq)
    
#     # Write the ccomverted sequences
#     write_seq_dict_to_file(sequence_dict, output_file)