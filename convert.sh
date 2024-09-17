#!/bin/bash

#   convert.sh
#   Version 1.0 10th September 2024

#   This shell script allows an Open Babel conversion to be carried out on the Python Flask app from the command line.

# Required positional arguments:
# arg1 is the file name, not including the extension (must not contain spaces, e.g., beta-carotene), or the word batch if multiple files are to be converted
# arg2 is the input (convert from) file extension (e.g., smi); for batch runs, files with other extensions are ignored
# arg3 is the output (convert to) file extension (e.g., cml)
# arg4 is the absolute path to the directory containing the input file or files (e.g., /Users/username/files/)
# arg5 is the absolute path to the download location (e.g., /Users/username/files/output/)

# Optional arguments must be placed before the required arguments:
# -h displays help (other arguments are ignored)
# -i value is a string of input (read) Open Babel option flags (e.g., a or aS, where a and S are individual options)
# -l Boolean value is true if a log file is required, false otherwise (default: true)
# -o value is a string of output (write) Open Babel option flags (e.g., m or 3m)
# -g values are Gen2D or Gen3D, for Open Babel calculation of 2D or 3D atomic coordinates respectively (ignored if not appropriate to the file format)
# -s values are fastest, fast, medium, better or best, and are used in conjuncion with -g (default: medium)

# To run in a terminal window, assuming that the current directory contains this file:
# sh convert.sh [-h] [-l <Boolean>] [-i <string>] [-o <string>] [-g <string>] [-s <string>] arg1 arg2 arg3 arg4 arg5

# For example:
# sh convert.sh -l false -i aS -o 3m beta-carotene smi cml /Users/username/files/ /Users/username/files/output/

from_flags=""
to_flags=""
help=false
log=true
gen="neither"
gen_opt="medium"

# Retrieve optional arguments
while getopts "hi:l:o:g:s:" flag
  do
    case $flag in
      h) help=true
         echo "
Required positional arguments:
  arg1 is the file name, not including the extension (must not contain spaces, e.g., beta-carotene), or the word batch if multiple files are to be converted
  arg2 is the input (convert from) file extension (e.g., smi); for batch runs, files with other extensions are ignored
  arg3 is the output (convert to) file extension (e.g., cml)
  arg4 is the absolute path to the directory containing the input file or files (e.g., /Users/username/files/)
  arg5 is the absolute path to the download location (e.g., /Users/username/files/output/)

Optional arguments must be placed before the required arguments:
  -h displays help (other arguments are ignored)
  -i is a string of input (read) Open Babel option flags (e.g., a or aS, where a and S are individual options)
  -l Boolean value is true if a log file is required, false otherwise (default: true)
  -o is a string of output (write) Open Babel option flags (e.g., m or 3m)
  -g values are Gen2D or Gen3D, for Open Babel calculation of 2D or 3D atomic coordinates respectively (ignored if not appropriate to the file format)
  -s values are fastest, fast, medium, better or best, and are used in conjuncion with -g (default: medium)

To run in a terminal window, assuming that the current directory contains this file:
  sh convert.sh [-h] [-l <Boolean>] [-i <string>] [-o <string>] [-g <string>] [-s <string>] arg1 arg2 arg3 arg4 arg5

For example:
  sh convert.sh -l false -i aS -o 3m beta-carotene smi cml /Users/username/files/ /Users/username/files/output/
";;
      i) from_flags=$OPTARG;;
      l) log=$OPTARG;;
      o) to_flags=$OPTARG;;
      g) gen=$OPTARG;;
      s) gen_opt=$OPTARG;;
    esac
  done

convert() {
  # Upload input file, convert and create log file
  curl -F 'from='"$arg2" -F 'to='"$arg3" -F 'from_flags='"$from_flags" -F 'to_flags='"$to_flags" -F 'coordinates='"$gen" -F 'coordOption='"$gen_opt" -F 'file=@'"$arg4$1"'.'"$arg2" http://127.0.0.1:5000/conv/
  #curl -F 'from='"$arg2" -F 'to='"$arg3" -F 'from_flags='"$from_flags" -F 'to_flags='"$to_flags" -F 'coordinates='"$gen" -F 'coordOption='"$gen_opt" -F 'file=@'"$arg4$1"'.'"$arg2" https://psdidev2.azurewebsites.net/conv/

  # Download output file
  curl -o $arg5$1.$arg3 http://127.0.0.1:5000/static/downloads/$1.$arg3
  #curl -o $arg5$1.$arg3 https://psdidev2.azurewebsites.net/static/downloads/$1.$arg3

  if [ $log = true ]; then
    # Download log file
    curl -o $arg5$1.log.txt http://127.0.0.1:5000/static/downloads/$1.log.txt
    #curl -o $arg5$1.log.txt https://psdidev2.azurewebsites.net/static/downloads/$1.log.txt
  fi

  # Delete input, output and log files from the server
  curl -d "filepath=static/uploads/$1.$arg2" http://127.0.0.1:5000/del/
  curl -d "filepath=static/downloads/$1.$arg3" http://127.0.0.1:5000/del/
  curl -d "filepath=static/downloads/$1.log.txt" http://127.0.0.1:5000/del/
  #curl -d "filepath=static/uploads/$1.$arg2" https://psdidev2.azurewebsites.net/del/
  #curl -d "filepath=static/downloads/$1.$arg3" https://psdidev2.azurewebsites.net/del/
  #curl -d "filepath=static/downloads/$1.log.txt" https://psdidev2.azurewebsites.net/del/;
}

if [ $help = false ]; then
  # Retrieve required arguments
  arg1=${@:OPTIND:1}
  arg2=${@:OPTIND+1:1}
  arg3=${@:OPTIND+2:1}
  arg4=${@:OPTIND+3:1}
  arg5=${@:OPTIND+4:1}

  if [ $arg1 = "batch" ]; then
    # Multiple files
    for file in "$arg4"*; do
      if [ -f "$file" ]; then
        filename=$(basename -- "$file")
        extension="${filename##*.}"
        filename="${filename%.*}"

        if [ $extension = $arg2 ]; then
          convert $filename
        fi
      fi
    done
  else
    # Single file
    convert $arg1
  fi
fi

