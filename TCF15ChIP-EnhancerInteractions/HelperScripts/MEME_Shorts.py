import subprocess
from pybedtools import BedTool

"""Collection of slightly easier access to the MEME suite tools."""


def streme(foreground, background='', sequence='', file_tag='', out_dir='', nmotifs=None, streme_src='streme'):
    """Run MEME's STREME to discover motifs that are enriched in the foreground compared to background."""
    if type(foreground) == str:
        foreground = BedTool(foreground)
    if background and type(background) == str:
        background = BedTool(background)
        background_path = out_dir + file_tag + "_BG.fa"
        background.sequence(fi=sequence, fo=background_path, name=False)
        background_flag = ' --n ' + background_path
    else:
        background_flag = ''
    if nmotifs:
        nmotifs_flag = ' --nmotifs ' + str(nmotifs)
    else:
        nmotifs_flag = ''
    foreground_path = out_dir + file_tag + "_FG.fa"
    foreground.sequence(fi=sequence, fo=foreground_path, name=False)
    streme_cmd = streme_src+" --p " + foreground_path + background_flag + " --o " + out_dir + file_tag + " --minw 6" + nmotifs_flag
    print(streme_cmd)
    subprocess.call(streme_cmd, shell=True)

