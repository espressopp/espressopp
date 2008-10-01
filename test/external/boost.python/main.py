# append the local path to be able to start in the build directory
import sys
sys.path.append(".")

import hello

hello.HelloWorld().printMessage()
