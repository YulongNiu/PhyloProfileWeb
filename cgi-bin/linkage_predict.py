#!/usr/local/bin/python

print "Content-Type: text/html\n"

# tmpData path
tempPath = '/var/www/html/phyloprofile/tmpData/'

# Debugging, and should be removed after publish the webpage.
import cgitb
cgitb.enable()

import cgi, os
from set_uni_file import RandomName
from rpy2.robjects import r, globalenv, FloatVector, StrVector
from rpy2.robjects.packages import importr
importr('hwriter')
importr('PhyloProfile')


######################Python_Funciton###############
def GetRFilePath(folderName, fileName):
    """
    'GetRFilePath' is used to generate a R file path in python CGI
    """
    filepwd = StrVector(folderName + '/'+ fileName)
    filepwd = r['paste'](filepwd, collapse = '')
    return filepwd
####################################################

fnDir = RandomName('phylopred')
fn = tempPath + fnDir
if os.path.exists(fn):
    fnDir = RandomName('phylopred')
    fn = tempPath + fnDir
os.mkdir(fn)

# return message whethe the upload data is correct
message = []
# get uploaded file
form = cgi.FieldStorage()

# A nested FieldStorage instance holds the file
fileitem = form['candidateList']
if fileitem.filename:
    # strip leading path from file name to avoid directory traversal attacks
    filen = os.path.basename(fileitem.filename)
    open(fn + '/' + filen, 'wb').write(fileitem.file.read())
else:
    message.append('No file was uploaded.\n')



#########################Process_data##################
if len(message) == 0:
    # source the R code. After that I just need to load the package.
    r.source('phylo_getlinkages.R')
    r.load('top400Viz.RData')
    r.load('profileViz.RData')
    top400Viz = r['top400Viz']
    wholePhyloDataNet = r['wholePhyloDataNet']
    kingdomCol = r['kingdomCol']
    
    # read in file
    if filen.split('.')[-1] == 'csv':
        # for csv file
        batArgu = r['read.csv'](fn + '/' + filen, stringsAsFactors = False)
        geneList = batArgu.rx(True, 1)
        geneColVec = batArgu.rx(True, 2)
    elif filen.split('.')[-1] == 'txt':
        # for txt file
        batArgu = r['read.table'](fn + '/' + filen, sep = '\t', header = True, stringsAsFactors = False, **{'comment.char': ''})
        geneList = batArgu.rx(True, 1)
        geneColVec = batArgu.rx(True, 2)
    else:
        message.append('Only "txt" or "csv" file format is allowed.\n')

        
    # No warning message
    if len(message) == 0:

        ##~~~~~~~~~~~~~~~~~~~~~~~~plot phyloprofile~~~~~~~~~~~~~~~~~~~
        # select profiles
        profileMat = r['GetProfile'](geneList, wholePhyloDataNet)
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
        globalenv['linksMat'] = linksMat
        linksMat = r("hwrite(linksMat, center = TRUE, row.bgcolor = '#a7c942', col.bgcolor = list('From' = '#a7c942', 'To' = '#a7c942'), row.style = list('font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'), col.style = list('From' = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em', 'To' = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'), br = TRUE, table.class = 'table')")
        linksMat = tuple(linksMat)[0]
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ##~~~~~~~~~~~~~~~~~~~~~~Generate_HTML_CSS_file~~~~~~~~~~~~~~~~~~~~~
        # read html template
        htmltemp = open('/var/www/cgi-bin/' + 'phylo_linkages.html').read()
        htmlReturn = htmltemp %(profileFigObj, cormatrixFigObj, linksMat, fnDir)
        # beware of the path!!!!!!!!!!!!!
        # write html index
        f = open(fn + '/index.html', 'w')
        f.write(htmlReturn + '\n')
        f.close()

        # # read css template
        # csstemp = open('/var/www/html/metaker/CSS/style2.css').read()
        # f = open(fn+'/style2.css', 'w')
        # f.write(csstemp + '\n')
        # f.close()

        # # copy images
        # imgdir = '/var/www/html/metaker/Pic/'
        # os.system('cp ' + imgdir + 'LOGO_PIC.png ' + fn)
        # os.system('cp ' + imgdir + 'LOGO_METAKER.png ' + fn)
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ##~~~~~~~~~~~~~~~~~~~~~~~~tar_Results~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        tarcom = 'tar -zcvf ' + fn + '.tar.gz ' '-C '+ tempPath + ' ' + fnDir
        os.system(tarcom + '>/dev/null')
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        #####################Return_HMTL_file###################
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

    else:
        message = StrVector(message)
        message = r.hwrite(message)
        message = tuple(message)[0]
        print """\
        <html><body>
        %s
        </body></html>
        """ % (message, )
else:
    message = StrVector(message)
    message = r.hwrite(message)
    message = tuple(message)[0]
    print """\
    <html><body>
    %s
    </body></html>
    """ % (message, )
