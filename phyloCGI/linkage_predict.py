#!/usr/local/bin/python

print("Content-Type: text/html\n")

# tmpData path
tempPath = '/var/www/html/phyloprofile/tmpData/'

# Initial warning message
message = ''

# Debugging, and should be removed after publish the webpage.
import cgitb
cgitb.enable()

import cgi, os, sys, math
from set_uni_file import RandomName
from rpy2.robjects import r
from rpy2.robjects.vectors import StrVector, IntVector, FloatVector
from rpy2.robjects.packages import importr

######################Python_Funciton###############
def GetRFilePath(folderName, fileName):
    """
    'GetRFilePath' is used to generate a R file path in python CGI
    """
    filepwd = StrVector(folderName + fileName)
    filepwd = r['paste'](filepwd, collapse = '')
    return filepwd
####################################################


#~~~~~~~~~~~ check threshold and input file~~~~~~~~
form = cgi.FieldStorage()

# get top number threshold
topNum = form.getfirst('topNum', '')
topNum = math.ceil(float(topNum))
topNum = int(topNum)

if (not 1 <= topNum <= 500):
    message += 'The threshold should be set between 1 and 500.\n'
    print('<html><body>%s</body></html>' % message)
    sys.exit(0)

# get phylogenetic plot para
phyloGeneNameSize = float(form.getfirst('phyloGeneNameSize', ''))
# get correlation plot para
corGeneNameSize = float(form.getfirst('corGeneNameSize', ''))

# A nested FieldStorage instance holds the file
fileitem = form['candidateList']
if fileitem.filename:
    # strip leading path from file name to avoid directory traversal attacks
    filen = fileitem.filename
else:
    message += 'No file was uploaded.\n'
    print('<html><body>%s</body></html>' % message)
    sys.exit(0)

fileSuffix = filen.split('.')[-1]

if fileSuffix not in ['csv', 'txt']:
    message += 'Only "txt" or "csv" file format is allowed.\n'
    print('<html><body>%s</body></html>' % message)
    sys.exit(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~get file~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~build tmp folder name~~~~~~~~~
while True:
    fnDir = RandomName('phylopred')
    fn = tempPath + fnDir + '/'
    if os.path.exists(fn) is not True:
        # get an unique name 
        os.mkdir(fn)
        open(fn + filen, 'wb').write(fileitem.file.read())
        break
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read in file
if fileSuffix == 'csv':
    # for csv file
    batArgu = r['read.csv'](fn + filen, stringsAsFactors = False)
else:
    # for txt file
    batArgu = r['read.table'](fn + filen, sep = '\t', header = True,
                              stringsAsFactors = False, **{"comment.char": ""})

# check the number of input genes should be at least two
if len(batArgu.rx(True, 1)) < 2:
    message += 'The input gene number should be at least two.\n'
    print('<html><body>%s</body></html>' % message)
    sys.exit(0)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#########################Process_data##################
# load library
importr('hwriter')
importr('PhyloProfile')
importr('PhyloProfileSuppl')

# I may need to use a database instead of load RData
r['load']('top500List.RData')
r['load']('wholePhyloDataAnno.RData')
r['load']('wholeProfile.RData')
r['load']('kingdomCol.RData')
r['load']('geneAnno.RData')
r['load']('phyloSpe.RData')
r['load']('kingdomAnno.RData')
top500List = r['top500List']
wholePhyloDataAnno = r['wholePhyloDataAnno']
wholeProfile = r['wholeProfile']
kingdomCol = r['kingdomCol']
geneAnno = r['geneAnno']
phyloSpe = r['phyloSpe']
kingdomAnno = r['kingdomAnno']

# retrieve vec
geneList = batArgu.rx(True, 1)
geneColVec = batArgu.rx(True, 2)
linkColVec = batArgu.rx(True, 3)

##~~~~~~~~~~~~~~~~~~~~~~~~~retrieve profile data~~~~~~~~~~~
# select profiles
profileMat = r['GetProfile'](geneList, wholeProfile)

# annotation
if batArgu.ncol == 4:
    # transfer gene anno to rownames
    usrGeneName = batArgu.rx(True, 4)
    geneMatchIdx = r['match'](profileMat.rownames, geneList)
    profileMat.rownames = usrGeneName.rx(geneMatchIdx)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~plot phyloprofile~~~~~~~~~~~~~~~~~~~
# set names of gene colors vector
geneColVec.names = geneList

# set profile plot path
profileFigPdfPath = GetRFilePath(fn, 'profilePlot.pdf')
profileFigJpgPath = GetRFilePath(fn, 'profilePlot.jpg')
        
r['pdf'](profileFigPdfPath, width = 10)
profileFig = r['PlotPhyloProfile'](profileMat,
                                   geneNameSize = phyloGeneNameSize,
                                   speCol = kingdomCol,
                                   geneCol = geneColVec,
                                   classCol = kingdomAnno,
                                   widthsShinkage = FloatVector([0.7, 0.9, 0.3, 7, 2]),
                                   **{"legend.position": "left"})

r['dev.off']()

os.system('convert -density 100 ' + ''.join(list(profileFigPdfPath)) +' ' + ''.join(list(profileFigJpgPath)) + ' >/dev/null')

profileFigObj = r("hwriteImage('profilePlot.jpg', center = TRUE)")
profileFigObj = tuple(profileFigObj)[0]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~plot cormatrix~~~~~~~~~~~~~~~~~~~
# set profile plot path
cormatrixFigPdfPath = GetRFilePath(fn, 'cormatrixPlot.pdf')
cormatrixFigJpgPath = GetRFilePath(fn, 'cormatrixPlot.jpg')

r['pdf'](cormatrixFigPdfPath)
cormatrixFig = r['PlotPhyloCor'](profileMat,
                                 geneNameSize = corGeneNameSize,
                                 geneCol = geneColVec,
                                 widthsShinkage = FloatVector([0.9, 0.9, 0.3, 7]),
                                 showCorVal = False)
r['dev.off']()

os.system('convert -density 100 ' + ''.join(list(cormatrixFigPdfPath)) + ' ' + ''.join(list(cormatrixFigJpgPath)) + ' >/dev/null')

cormatrixFigObj = r("hwriteImage('cormatrixPlot.jpg', center = TRUE)")
cormatrixFigObj = tuple(cormatrixFigObj)[0]

# write correlation matrix
corMat = r['GetPhyloCorMat'](profileMat)
corMatpwd = GetRFilePath(fn, 'correlation_matrix.csv')
r['write.csv'](corMat, corMatpwd)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~interaction matrix~~~~~~~~~~~~~~~~~~~~
linksMat = r['GetTopLink'](geneIDs = geneList,
                             linkData = top500List,
                             annoVec = wholePhyloDataAnno,
                             threshold = topNum)
linksMatpwd = GetRFilePath(fn, 'predicted_linakges.csv')
r['write.csv'](linksMat, linksMatpwd)
linksMatObj = r['hwrite'](linksMat, center = True, br = True,
                          **{"row.bgcolor": "#a7c942",
                             "col.bgcolor": r['list'](From = '#a7c942', To = '#a7c942'),
                             "row.style": r['list']('font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'),
                             "col.style": r['list'](From = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em', To = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'),
                             "table.class": "'table'"})
linksMatObj = tuple(linksMatObj)[0]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~circos plot~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check link colours
checkColList = r['CheckLinkCol'](geneList, linkColVec, geneAnno)
geneList = checkColList.rx2('checkGeneVec')
linkColVec = checkColList.rx2('checkLinkCol')
geneSym = checkColList.rx2('checkGeneSym')
wm = checkColList.rx2('wm')

if len(wm) == 1:
    circosFigObj = tuple(wm)[0]
elif len(wm) == 0 and len(geneList) > 7:
    circosFigObj = 'The number of candidate genes for Circos plot should be no more than 7.\n'
elif len(wm) == 0 and len(geneList) <= 7:
    # copy and compress circos folder
    os.system('cp circosConfig.tar.gz ' + fn + ' >/dev/null')
    os.system('tar -zxvf ' + fn + 'circosConfig.tar.gz -C ' + fn + ' >/dev/null')
    os.system('rm ' + fn + 'circosConfig.tar.gz' + ' >/dev/null')

    # generate circos files
    ftMat = linksMat.rx(True, IntVector((1, 3, 5)))
    r['writeCircos'](geneList, ftMat, geneAnno, phyloSpe, wholeProfile, savePath = fn + 'circosConfig/phylo/')

    # generate circos config
    r['writeConf']('phylo/', geneList, geneSym, linkColVec, fn + 'circosConfig/')

    # circos plot
    os.system('circos-0.67-7/bin/circos -conf ' + fn + 'circosConfig/circosConf.conf ' + '-outputdir ' + fn + ' -outputfile ' + 'circosPlot' + ' >/dev/null')
    os.system('convert -resize 700x700 ' + fn + 'circosPlot.png ' + fn + 'circosPlotWeb.png' + ' >/dev/null')
    circosFigObj = r['hwriteImage']('circosPlotWeb.png', center = True)
    circosFigObj = tuple(circosFigObj)[0]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d3 network~~~~~~~~~~~~~~~~~~~~~~~~~
if len(wm) == 1:
    d3FigObj = tuple(wm)[0]
elif len(wm) == 0 and len(geneList) > 7:
    d3FigObj = 'The number of candidate genes for D3network plot should be no more than 7.\n'
elif len(wm) == 0 and len(geneList) <= 7:
    # selection and annotation ftmat
    ftMat = linksMat.rx(True, IntVector((1, 3, 5)))
    annoftMat = r['Annoft'](geneList, ftMat, geneAnno)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~Generate_HTML_CSS_file~~~~~~~~~~~~~~~~~~~~~
# read html template
htmltemp = open('/var/www/cgi-bin/phyloCGI/' + 'phylo_linkages.html').read()
htmlReturn = htmltemp %(topNum, profileFigObj, cormatrixFigObj, circosFigObj, linksMatObj, fnDir)
# beware of the path!!!!!!!!!!!!!
# write html index
f = open(fn + 'index.html', 'w')
f.write(htmlReturn + '\n')
f.close()

# # read css template
# csstemp = open('/var/www/html/metaker/CSS/style2.css').read()
# f = open(fn + 'style2.css', 'w')
# f.write(csstemp + '\n')
# f.close()

# # copy images
# imgdir = '/var/www/html/metaker/Pic/'
# os.system('cp ' + imgdir + 'LOGO_PIC.png ' + fn)
# os.system('cp ' + imgdir + 'LOGO_METAKER.png ' + fn)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~tar_Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tarcom = 'tar -zcvf ' + fn[ :-1] + '.tar.gz ' '-C '+ tempPath + ' ' + fnDir
os.system(tarcom + ' >/dev/null')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~Return_HMTL_file~~~~~~~~~~~~~~~~~~~~~~~
# print htmlReturn
print("""\
<html>
<head>
<meta http-equiv="Refresh" content="2;url=/phyloprofile/tmpData/%s/index.html">
<base target="_blank"/>
</head>
<body>
<p>
New address in 2s.
</p>
</body>
</html>
""") % fnDir
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
