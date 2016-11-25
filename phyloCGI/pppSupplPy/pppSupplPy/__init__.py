"""
'RandomName' is used to randomly generate a name.
The length of file name will be set to 20.
"""
import random, os
def RandomName(fileType):
    # low 26 letters
    letter = [chr(i) for i in range(97, 123)]
    # up 26 letters
    LETTER = [chr(i) for i in range(65, 91)]
    # 0-9
    num = list(range(10))

    # randomly choose name
    name = [str(random.choice(letter + LETTER + num)) for i in range(20)]
    fullName = fileType + '_' + ''.join(name)

    return fullName


"""
'GetRFilePath' is used to generate a R file path in python CGI
"""
from rpy2.robjects import r
from rpy2.robjects.vectors import StrVector
def GetRFilePath(folderName, fileName):
    """
    'GetRFilePath' is used to generate a R file path in python CGI
    """
    filepwd = StrVector(folderName + fileName)
    filepwd = r['paste'](filepwd, collapse = '')
    return filepwd
