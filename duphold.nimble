import ospaths
template thisModuleFile: string = instantiationInfo(fullPaths = true).filename

when fileExists(thisModuleFile.parentDir / "src/duphold.nim"):
  # In the git repository the Nimble sources are in a ``src`` directory.
  import src/dupholdpkg/version as _
else:
  # When the package is installed, the ``src`` directory disappears.
  import dupholdpkg/version as _

# Package

version       = dupholdVersion
author        = "Brent Pedersen"
description   = "find depth support for DUP/DEL/CNV calls that use PE/SR"
license       = "MIT"


# Dependencies

requires "docopt#0abba63", "genoiser >= 0.2.2", "hts >= 0.2.5"
srcDir = "src"
installExt = @["nim"]

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

