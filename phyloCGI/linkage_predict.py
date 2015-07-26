#!/usr/local/bin/python

print "Content-Type: text/html\n"

# tmpData path
tempPath = '/var/www/html/phyloprofile/tmpData/'

# Debugging, and should be removed after publish the webpage.
import cgitb
cgitb.enable()

import cgi, os
from set_uni_file import RandomName
from rpy2.robjects import r, StrVector, IntVector
from rpy2.robjects.packages import importr
importr('hwriter')
importr('PhyloProfile')
importr('PhyloProfileSuppl')


######################Python_Funciton###############
def GetRFilePath(folderName, fileName):
    """
    'GetRFilePath' is used to generate a R file path in python CGI
    """
    filepwd = StrVector(folderName + fileName)
    filepwd = r['paste'](filepwd, collapse = '')
    return filepwd
####################################################

fnDir = RandomName('phylopred')
fn = tempPath + fnDir + '/'
if os.path.exists(fn):
    fnDir = RandomName('phylopred')
    fn = tempPath + fnDir + '/'
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
    open(fn + filen, 'wb').write(fileitem.file.read())
else:
    message.append('No file was uploaded.\n')



#########################Process_data##################
if len(message) == 0:
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
    
    # read in file
    if filen.split('.')[-1] == 'csv':
        # for csv file
        batArgu = r['read.csv'](fn + filen, stringsAsFactors = False)
    elif filen.split('.')[-1] == 'txt':
        # for txt file
        batArgu = r['read.table'](fn + filen, sep = '\t', header = True, stringsAsFactors = False, **{"comment.char": ""})
    else:
        message.append('Only "txt" or "csv" file format is allowed.\n')

    # No warning message
    if len(message) == 0:

        # retrieve vec
        geneList = batArgu.rx(True, 1)
        geneColVec = batArgu.rx(True, 2)
        linkColVec = batArgu.rx(True, 3)

        # check row number should be more than 1
        if list(r['length'](geneList))[0] < 2:
            message.append('The input gene number should be at least two.\n')

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
        linksMatObj = r['hwrite'](linksMat, center = True, br = True, **{"row.bgcolor": "#a7c942", "col.bgcolor": r['list'](From = '#a7c942', To = '#a7c942'), "row.style": r['list']('font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'), "col.style": r['list'](From = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em', To = 'font-weight:bold; text-align:center; color:#413b62; padding-top:5px; padding-bottom:4px; font-size:1.1em'), "table.class": "'table'"})
        linksMatObj = tuple(linksMatObj)[0]
        ##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ##~~~~~~~~~~~~~~~~~~~~~~~~~circos plot~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # check link colours
        checkColList = r['CheckLinkCol'](geneList, linkColVec, geneAnno)
        wm = checkColList.rx2(4)
        if list(r['length'](wm))[0] == 0:
            geneList = checkColList.rx2(1)
            linkColVec = checkColList.rx2(2)
            geneSym = checkColList.rx2(3)

            # copy and compress circos folder
            os.system('cp circosConfig.tar.gz ' + fn + ' >/dev/null')
            os.system('tar -zxvf ' + fn + 'circosConfig.tar.gz -C ' + fn + ' >/dev/null')

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
