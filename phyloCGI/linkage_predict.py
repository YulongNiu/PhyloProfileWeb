#!/usr/local/bin/python

print "Content-Type: text/html\n"

# tmpData path
tempPath = '/var/www/html/phyloprofile/tmpData/'

# Initial warning message
message = ''

# Debugging, and should be removed after publish the webpage.
import cgitb
cgitb.enable()

import cgi, os, sys
from set_uni_file import RandomName
from rpy2.robjects import r
from rpy2.robjects.vectors import StrVector, IntVector
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


#~~~~~~~~~~~ check uploaded file~~~~~~~~
form = cgi.FieldStorage()

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
r['load']('top400Viz.RData')
r['load']('wholeProfile.RData')
r['load']('kingdomCol.RData')
r['load']('geneAnno.RData')
r['load']('phyloSpe.RData')
top400Viz = r['top400Viz']
wholeProfile = r['wholeProfile']
kingdomCol = r['kingdomCol']
geneAnno = r['geneAnno']
phyloSpe = r['phyloSpe']

# retrieve vec
geneList = batArgu.rx(True, 1)
geneColVec = batArgu.rx(True, 2)
linkColVec = batArgu.rx(True, 3)

##~~~~~~~~~~~~~~~~~~~~~~~~plot phyloprofile~~~~~~~~~~~~~~~~~~~
# select profiles
profileMat = r['GetProfile'](geneList, wholeProfile)
# set names of gene colors vector
geneColVec.names = geneList

# set profile plot path
profileFigPdfPath = GetRFilePath(fn, 'profilePlot.pdf')
profileFigJpgPath = GetRFilePath(fn, 'profilePlot.jpg')
        
r['pdf'](profileFigPdfPath)
profileFig = r['PlotPhyloProfile'](profileMat, speCol = kingdomCol, geneCol = geneColVec)
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
cormatrixFig = r['PlotPhyloCor'](profileMat, geneCol = geneColVec)
r['dev.off']()

os.system('convert -density 100 ' + ''.join(list(cormatrixFigPdfPath)) + ' ' + ''.join(list(cormatrixFigJpgPath)) + ' >/dev/null')

cormatrixFigObj = r("hwriteImage('cormatrixPlot.jpg', center = TRUE)")
cormatrixFigObj = tuple(cormatrixFigObj)[0]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
##~~~~~~~~~~~~~~~~~~~~~~~~~interaction matrix~~~~~~~~~~~~~~~~~~~~
linksMat = r['GetLinkages'](geneList, top400Viz)
linksMatpwd = GetRFilePath(fn, 'predicted_linakges.csv')
r['write.csv'](linksMat, linksMatpwd)
linksMatObj = r['hwrite'](linksMat, center = True, br = True,
                          **{"row.bgcolor": "#a7c942",
                             "col.bgcolor": r['list'](From = '#a7c942', To = '#a7c942'), "row.style": r['list']('font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'),
                             "col.style": r['list'](From = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em', To = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'),
                             "table.class": "'table'"})
linksMatObj = tuple(linksMatObj)[0]
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~~~~circos plot~~~~~~~~~~~~~~~~~~~~~~~~~~~
# check link colours
checkColList = r['CheckLinkCol'](geneList, linkColVec, geneAnno)
wm = checkColList.rx2(4)
if len(wm) == 0:
    geneList = checkColList.rx2(1)
    linkColVec = checkColList.rx2(2)
    geneSym = checkColList.rx2(3)
    
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
else:
    wm = r.hwrite(wm)
    circosFigObj = tuple(wm)[0]     
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##~~~~~~~~~~~~~~~~~~~~~~Generate_HTML_CSS_file~~~~~~~~~~~~~~~~~~~~~
# read html template
htmltemp = open('/var/www/cgi-bin/phyloCGI/' + 'phylo_linkages.html').read()
htmlReturn = htmltemp %(profileFigObj, cormatrixFigObj, circosFigObj, linksMatObj, fnDir)
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
print """\
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
""" % fnDir
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##################################################################
