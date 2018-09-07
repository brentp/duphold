# Package

version       = "0.0.1"
author        = "Brent Pedersen"
description   = "find depth support for DUP/DEL/CNV calls that use PE/SR"
license       = "MIT"


# Dependencies

requires "hts >= 0.2.0", "docopt", "genoiser"
srcDir = "src"

bin = @["duphold"]


skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "nim c --lineDir:on --debuginfo -r --threads:on tests/all"

#before test:
#  exec "c2nim src/hts/private/hts_concat.h"

task docs, "Builds documentation":
  mkDir("docs"/"duphold")
  #exec "nim doc2 --verbosity:0 --hints:off -o:docs/index.html  src/hts.nim"
  for file in listfiles("src/duphold"):
    if file.endswith("value.nim"): continue
    if splitfile(file).ext == ".nim":
      exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ file.changefileext("html").split("/", 1)[1] & " " & file

