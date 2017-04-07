#!/bin/sh -v
# $Id: gen_html.sh,v 1.1 1999/09/14 09:51:25 aimar Exp $
# Massages the HTML generated from Save as HTML
#
#

\cp ./sdlt-master.css html/DevGuide.css

\rm  html/welcome.html
\rm  html/index.html
\rm  html/contents.html

\ln -s DevGuide.html    html/welcome.html
\ln -s DevGuide.html    html/index.html
\ln -s DevGuide-1.html  html/contents.html
