#!/usr/bin/python

from paul.browser.browserwindow import *

testfile = "/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw"

app = QApplication (sys.argv)
w = BrowserWindow()
app.setMainWidget (w)

w.show()
sys.exit(app.exec_loop())
