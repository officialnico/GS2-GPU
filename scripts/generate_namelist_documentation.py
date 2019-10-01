#!/usr/bin/env python3

#Example usage from GK_HEAD_DIR
#> scripts/generate_namelist_documentation.py $(find ./src ./externals -name '*.[fF]??')
#Will produce namelist documentation for use with ford/bitbucket wiki
#Output can be found in docs_path+pages_path+namelists_path

from config_generation import constructFromFile
from get_all_matching_files import get_all_matching_files
docs_path = "docs/"
pages_path = "pages/"
namelists_path = "namelists/"
# If the following is set to true then we also generate
# a separate page per namelist along with the overview page.
# This can significantly clutter up the side bar so default
# to False here
one_page_per_namelist = False

if __name__ == "__main__":
    import sys
    import os

    # Handle case where no files passed as arguments
    if len(sys.argv) < 2:
        probablyFiles = get_all_matching_files(paths=["src/","externals/utils/"],
                                               patterns=["*.[fF]??"])
    else:
        probablyFiles = sys.argv[1:]

    # Check target directory exists
    target_dir = docs_path+pages_path+namelists_path

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    with open(docs_path+pages_path+namelists_path+"index.md", 'w') as ind:
        ind.write("title: Namelists\n")
        ind.write("\n")
        ind.write("\n[TOC]\n\n\n")        
        ind.write('<html><div class="divNamelistTable"></html>')
        ind.write("\n")        
        for fil in probablyFiles:
            tmp = constructFromFile(fileIn=fil)
            if tmp is None: continue
            for k, v in tmp.items():
                ff = pages_path+namelists_path+"auto_generated_"+v.moduleName+"_"+v.namelistName
                fullf = docs_path+ff
                with open(fullf+".txt", 'w') as f:
                    f.writelines("\n".join(v.toMarkdownTable()))

                if one_page_per_namelist:
                    with open(fullf+".md", 'w') as f:
                        f.write("title: {namelist}\n".format(namelist=v.namelistName))
                        f.write("\n")
                        f.write("\n".join(v.comment))
                        f.write("\n\n")
                        f.write('<html><div class="divNamelistTable"></html>')
                        f.write("\n\n")
                        f.write("{!"+ff+".txt"+"!}\n")
                        f.write("\n")
                        f.write('</div>')
                    
                ind.write("{!"+ff+".txt"+"!}\n")
                ind.write("\n")
        ind.write('</div>')

