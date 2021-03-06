#! /usr/local/python3/bin/python3

print("Content-Type: text/html\n")

## tmpData path
tempPath = '/var/www/html/phyloprofile/tmpData/'

## Initial warning message
message = ''

## Debugging, and should be removed after publish the webpage.
import cgitb
cgitb.enable()

import cgi, os, sys, math
from pppSupplPy import RandomName, GetRFilePath
from rpy2.robjects import r
from rpy2.robjects.vectors import IntVector, FloatVector
from rpy2.robjects.packages import importr
############### check threshold and input file ################
form = cgi.FieldStorage()

## get top number threshold
topNum = form.getfirst('topNum', '')
topNum = math.ceil(float(topNum))
topNum = int(topNum)

if (not 1 <= topNum <= 500):
    message += 'The threshold should be set between 1 and 500.\n'
    print('<html><body>%s</body></html>' % message)
    sys.exit(0)

## get blast evalue threshold
evalue = form.getfirst('evalue', '')
evalueObj = format(float(evalue), '.4f')

## get organism
org = form.getfirst('org', '')

## get phylogenetic plot para
phyloGeneNameSize = float(form.getfirst('phyloGeneNameSize', ''))
## get correlation plot para
corGeneNameSize = float(form.getfirst('corGeneNameSize', ''))

## A nested FieldStorage instance holds the file
fileitem = form['candidateList']
if fileitem.filename:
    ## strip leading path from file name to avoid directory traversal attacks
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
#########################################################

#########################get file###########################
##~~~~~~~~~~~~~~~~build tmp folder name~~~~~~~~~
while True:
    fnDir = RandomName('ppp')
    fn = tempPath + fnDir + '/'
    if os.path.exists(fn) is not True:
        # get an unique name
        os.mkdir(fn)
        open(fn + filen, 'wb').write(fileitem.file.read())
        break
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## read in file
if fileSuffix == 'csv':
    ## for csv file
    batArgu = r['read.csv'](fn + filen, stringsAsFactors = False)
else:
    ## for txt file
    batArgu = r['read.table'](fn + filen, sep = '\t', header = True,
                              stringsAsFactors = False, **{"comment.char": ""})

## check the number of input genes should be at least two
if len(batArgu.rx(True, 1)) < 2:
    message += 'The input gene number should be at least two.\n'
    print('<html><body>%s</body></html>' % message)
    sys.exit(0)
##########################################################

#########################Process_data##################
## load library
importr('PhyloProfile')
importr('pppSupplR')

## I may need to use a database instead of load RData
orgPrefix = org + '/' + org + '_'
r['load'](orgPrefix + 'top500List' + evalue + '.RData')
r['load'](orgPrefix + 'wholeGenomeAnno.RData')
r['load'](orgPrefix + 'wholeProfile' + evalue + '.RData')
r['load'](orgPrefix + 'geneAnno.RData')
r['load']('kingdomCol.RData')
r['load']('phyloSpe.RData')
r['load']('kingdomAnno.RData')
top500List = r['top500List']
wholeGenomeAnno = r['wholeGenomeAnno']
wholeProfile = r['wholeProfile']
kingdomCol = r['kingdomCol']
geneAnno = r['geneAnno']
phyloSpe = r['phyloSpe']
kingdomAnno = r['kingdomAnno']

## retrieve vec
geneList = batArgu.rx(True, 1)
geneColVec = batArgu.rx(True, 2)
linkColVec = batArgu.rx(True, 3)

##~~~~~~~~~~~~~~~~~~~~~~~~~retrieve profile data~~~~~~~~~~~
## select profiles
geneMatchIdx = r['match'](geneList,
                          wholeProfile.rownames,
                          nomatch = 0)
profileMat = wholeProfile.rx(geneMatchIdx, True)

## annotation
if batArgu.ncol == 4:
    ## transfer gene anno to rownames
    usrGeneName = batArgu.rx(geneMatchIdx.ro != 0, 4)
    profileMat.rownames = usrGeneName
    ## change gene colours names
    geneColVec.names = usrGeneName
else:
    geneColVec.names = geneList
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~plot phyloprofile~~~~~~~~~~~~~~~~~~~
## set profile plot path
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
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~plot cormatrix~~~~~~~~~~~~~~~~~~~
## set profile plot path
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

## write correlation matrix
corMat = r['GetPhyloCorMat'](profileMat)
corMatpwd = GetRFilePath(fn, 'correlation_matrix.csv')
r['write.csv'](corMat, corMatpwd)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~interaction matrix~~~~~~~~~~~~~~~~~~~~
linksMat = r['GetTopLink'](geneIDs = geneList,
                             linkData = top500List,
                             annoVec = wholeGenomeAnno,
                             threshold = topNum)
linksMatpwd = GetRFilePath(fn, 'predicted_linakges.csv')
r['write.csv'](linksMat, linksMatpwd, **{"row.names": False})
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~circosJS plot~~~~~~~~~~~~~~~~~~~~~~~~~~~
## check link colours
checkColList = r['CheckLinkCol'](geneList, linkColVec, geneAnno)
wm = checkColList.rx2('wm')

if len(wm) == 1:
    circosJSFigObj = tuple(wm)[0]
else:
    ## copy and compress circosJS folder
    os.system('cp circosConfig.tar.gz ' + fn + ' >/dev/null')
    os.system('tar -zxvf ' + fn + 'circosConfig.tar.gz -C ' + fn + ' >/dev/null')
    os.system('rm ' + fn + 'circosConfig.tar.gz' + ' >/dev/null')

    ## generate circosJS files
    ftMat = linksMat.rx(True, IntVector((1, 3, 5, 6)))
    r['writeCircosJS'](checkColList, ftMat, geneAnno, phyloSpe, wholeProfile, savePath = fn + 'circosConfig/')
    circosJSFigObj = ''
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d3 network~~~~~~~~~~~~~~~~~~~~~~~~~
if len(wm) == 1:
    d3FigObj = tuple(wm)[0]
else:
    ## selection and annotation ftmat
    annoftMat = r['read.csv'](fn + 'circosConfig/linkage.csv',
                              stringsAsFactor = False)
    annoftMat = annoftMat.rx(True, IntVector([4, 8, 9, 11]))

    ## transfer and plot de network
    d3ft = r['d3Transft'](annoftMat)
    d3Obj = r['d3PlotNet'](d3ft)

    ## write d3net
    r['writed3Net'](d3Obj, fileName = 'networkd3.html', savePath = fn)
    d3FigObj = r['d3ExtractNetEle'](fn + 'networkd3.html')

    ## remove old html file and js/css
    os.remove(fn + 'networkd3.html')
    os.system('rm -rf ' + fn + 'networkd3_files' + ' >/dev/null')

    d3FigObj = tuple(d3FigObj)[0]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~Generate_HTML_CSS_file~~~~~~~~~~~~~~~~~~~~~
## read html template
htmltemp = open('phylo_linkages.html').read()
replaceDic = {'topNum': topNum,
              'evalueObj': evalueObj,
              'circosJSFigObj': circosJSFigObj,
              'd3FigObj': d3FigObj,
              'fnDir': fnDir,
              'orgObj': org}
htmlReturn = htmltemp.format(**replaceDic)
## beware of the path!!!!!!!!!!!!!
## write html index
f = open(fn + 'index.html', 'w')
f.write(htmlReturn + '\n')
f.close()
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~tar_Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tarcom = 'tar -zcvf ' + fn[ :-1] + '.tar.gz ' + fn
os.system(tarcom + ' >/dev/null')
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~Return_HMTL_file~~~~~~~~~~~~~~~~~~~~~~~
## print htmlReturn
replaceDirDict = {'outdir':fnDir}
print("""\
<html>
<head>
<meta http-equiv="Refresh" content="2;url=/phyloprofile/tmpData/{outdir}/index.html">
<base target="_blank"/>
</head>
<body>
<p>
New address in 2s. Please refresh this page if no response in a long time.
</p>
</body>
</html>
""".format(**replaceDirDict))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
