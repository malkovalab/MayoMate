# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from mayo import csvDataFrameImport
from mayo.parental_reference_base import createParentalReferenceGenome
import mayo.settings.config as cfg

if __name__ == "__main__":

    original_ref_path = '/content/S288C_reference_sequencewoUra3wUra329wNAT.fa'
    new_ref_path =      '/content/S288C_reference_sequencewoUra3wUra329wNAT_mod2.fa'

    frames_list_parents = csvDataFrameImport(cfg.config["parents_path"], kind="parental")
    assert len(frames_list_parents) == 2, "There should be exactly 2 parents in the parents file"
    
    parent1, parent2 = frames_list_parents[0:2]

    createParentalReferenceGenome(parent1.df, parent2.df, original_ref_path, new_ref_path)