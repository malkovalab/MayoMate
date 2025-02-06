# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com

from .TableSNP import TableSNP
import logging

logger = logging.getLogger('__main__.' + __name__)

#should probably rename to SNPsTableParents
class ParentalTable(TableSNP):
    pass