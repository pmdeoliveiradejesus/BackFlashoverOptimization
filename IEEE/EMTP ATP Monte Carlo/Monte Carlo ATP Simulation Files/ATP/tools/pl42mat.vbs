set x=createobject("wscript.shell")
x.run "runGTPPL.exe bfr.pl4"
wscript.sleep 2000
x.sendkeys"matlab all"
x.sendkeys"{enter}"
x.sendkeys"{enter}"
x.sendkeys"{enter}"
wscript.sleep 2000
x.sendkeys"stop"
x.sendkeys"{enter}"
x.sendkeys"{enter}"
x.sendkeys"{enter}"