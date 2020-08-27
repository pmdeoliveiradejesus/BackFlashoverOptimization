set x=createobject("wscript.shell")
x.run "runATP.exe BFR.atp" 
wscript.sleep 2000
x.sendkeys"{enter}"
