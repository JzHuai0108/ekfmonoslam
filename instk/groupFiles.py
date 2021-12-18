"""Move .m files in instk to corresponding subfolders as in the package downloaded from the instk website."""
import os
import shutil


def filesInFolder(mypath):
    f = []
    for (dirpath, dirnames, filenames) in os.walk(mypath):
        for fn in filenames:
            f.append(os.path.join(dirpath, fn))
        break
    return f


def findFileInFolder(mypath, fn):
    f = []
    for (dirpath, dirnames, filenames) in os.walk(mypath):
        if fn in filenames:
            f.append(os.path.join(dirpath, fn))
    return f


localfolder = '/docker/ekfmonoslam/instk'
remotefolder = '/Downloads/MRepo'

files = filesInFolder(localfolder)

for f in files:
    name = os.path.basename(f)
    matches = findFileInFolder(remotefolder, name)

    if len(matches) == 1:
        dest = matches[0].replace(remotefolder, localfolder, 1)
        print('will move {} to {}'.format(f, dest))
        destFolder = os.path.dirname(dest)
        pathexist = os.path.isdir(destFolder)
        if not pathexist:
            os.makedirs(destFolder, exist_ok=True)
            print('Path {} created.'.format(destFolder))
        shutil.move(f, dest)
    elif len(matches) == 0:
        print('No matches for {}.'.format(name))
    else:
        print("Multiple matches for {}:\n{}".format(name, matches))


