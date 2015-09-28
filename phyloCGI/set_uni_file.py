"""
'RandomName' is used to randomly generate a name.
The length of file name will be set to 20, first of which will be capital letters.
"""
import random, os
def RandomName(fileType):
    # low 26 letters
    letter = [chr(i) for i in range(97, 123)]
    # up 26 letters
    LETTER = [chr(i) for i in range(65, 91)]
    # 0-9
    num = range(10)

    # randomly choose head (must be letters) and body
    head = [random.choice(LETTER) for i in range(1)]
    body = [str(random.choice(letter + LETTER + num)) for i in range(19)]

    name = ''.join(head + body)
    fullName = fileType + '_' + name
    return(fullName)


"""
'GetRFilePath' is used to generate a R file path in python CGI
"""
def GetRFilePath(folderName, fileName):
    """
    'GetRFilePath' is used to generate a R file path in python CGI
    """
    filepwd = StrVector(folderName + fileName)
    filepwd = r['paste'](filepwd, collapse = '')
    return filepwd
