import os

file_to_split = "whole.sdf"
direc = "chunks"

# each file will have at least this many lines
lines_per_file = 500000

try:
    os.mkdir("chunks")
except FileExistsError:
    pass

sname, ext = file_to_split.split(".")

# https://stackoverflow.com/questions/16840554/reading-a-line-from-file-without-advancing-pythonic-approach
def peek_line(f):
    pos = f.tell()
    line = f.readline()
    f.seek(pos)
    return line

with open(file_to_split) as main_file:
    fcount = 0
    while True:
        fname = direc + "/" + sname + str(fcount) + "." + ext
        lcount = 1
        peekline = main_file.readline()
        if peekline == "":
            break

        with open(fname, "w") as new_file:
            new_file.write(peekline)
            for line in main_file:
                new_file.write(line)
                lcount += 1
                if lcount > 500000 and line[:4] == "$$$$":
                    break;
        print(f"{fname} was created")
        fcount += 1;






