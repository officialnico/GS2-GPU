#!/usr/bin/env python3
"""Shims and stubs needed to monkey-patch FORD to allow for linking
to the GS2 namelists pages

"""

from operator import attrgetter
import os

import ford
from config_generation import (
    constructFromFile,
    config as Namelist,
    var as NamelistVariable,
)


NAMELIST_PAGE_URL = "{}/page/namelists/index.html#{}"


def get_namelist_url(self):
    """Get the url for a namelist
    """
    return NAMELIST_PAGE_URL.format(self.base_url, self.namelistName)


def get_namelist_variable_url(self):
    """Get the url for a specific variable in a namelist
    """
    anchor = "{}-{}".format(self.namelist, self.name)
    return NAMELIST_PAGE_URL.format(anchor)


def interceptor_main(proj_data, proj_docs, md):
    """
    Main driver of FORD with GS2-specific patches
    """

    if proj_data["relative"]:
        proj_data["project_url"] = "."

    # Parse the files in your project
    project = ford.fortran_project.Project(proj_data)
    if len(project.files) < 1:
        raise RuntimeError(
            "No source files with appropriate extension found in specified directory."
        )

    # Convert the documentation from Markdown to HTML. Make sure to properly
    # handle LateX and metadata.
    if proj_data["relative"]:
        project.markdown(md, "..")
    else:
        project.markdown(md, proj_data["project_url"])

    project.correlate()

    # This block is GS2 specific. We need to tell FORD about our
    # namelists
    project.namelists = []
    for file in project.files:
        tmp = constructFromFile(file.path)
        if tmp is None:
            continue
        for namelist in tmp.values():
            key = namelist.namelistName
            namelist.variables = namelist.varList

            # Make sure key is unique
            while key in map(attrgetter("name"), project.namelists):
                key += "_"

            for variable in namelist.variables:
                variable.namelist = key
                variable.anchor = "{}-{}".format(key, variable.name)

            namelist.base_url = (
                ".." if proj_data["relative"] else proj_data["project_url"]
            )
            namelist.name = key
            project.namelists.append(namelist)

    if proj_data["relative"]:
        project.make_links("..")
    else:
        project.make_links(proj_data["project_url"])

    # Convert summaries and descriptions to HTML
    if proj_data["relative"]:
        ford.sourceform.set_base_url(".")
    if "summary" in proj_data:
        proj_data["summary"] = md.convert(proj_data["summary"])
        proj_data["summary"] = ford.utils.sub_links(
            ford.utils.sub_macros(
                ford.utils.sub_notes(proj_data["summary"]), proj_data["project_url"]
            ),
            project,
        )
    if "author_description" in proj_data:
        proj_data["author_description"] = md.convert(proj_data["author_description"])
        proj_data["author_description"] = ford.utils.sub_links(
            ford.utils.sub_macros(
                ford.utils.sub_notes(proj_data["author_description"]),
                proj_data["project_url"],
            ),
            project,
        )

    proj_docs_ = ford.utils.sub_links(
        ford.utils.sub_macros(
            ford.utils.sub_notes(proj_docs), proj_data["project_url"]
        ),
        project,
    )

    # Process any pages
    if "page_dir" in proj_data:
        page_tree = ford.pagetree.get_page_tree(
            os.path.normpath(proj_data["page_dir"]), md
        )
        print()
    else:
        page_tree = None

    proj_data["pages"] = page_tree

    # This block is GS2 specific. FORD also does this monkey-patching
    # itself, inside PagetreePage.render. Easier to just do this
    # ourselves here though
    if proj_data["relative"]:
        for namelist in project.namelists:
            namelist.base_url = "../.."

    # Produce the documentation using Jinja2. Output it to the desired location
    # and copy any files that are needed (CSS, JS, images, fonts, source files,
    # etc.)

    docs = ford.output.Documentation(proj_data, proj_docs_, project, page_tree)
    docs.writeout()
    print("")


def interceptor_sub_links(string, project):
    """
    Replace links to different parts of the program, formatted as
    [[name]] or [[name(object-type)]] with the appropriate URL. Can also
    link to an item's entry in another's page with the syntax
    [[parent-name:name]]. The object type can be placed in parentheses
    for either or both of these parts.

    Patched with the GS2-specific 'namelist' knowledge
    """
    LINK_TYPES = {
        "module": "modules",
        "type": "types",
        "procedure": "procedures",
        "subroutine": "procedures",
        "function": "procedures",
        "proc": "procedures",
        "file": "allfiles",
        "interface": "absinterfaces",
        "absinterface": "absinterfaces",
        "program": "programs",
        "block": "blockdata",
        "namelist": "namelists",
    }

    SUBLINK_TYPES = {
        "variable": "variables",
        "type": "types",
        "constructor": "constructor",
        "interface": "interfaces",
        "absinterface": "absinterfaces",
        "subroutine": "subroutines",
        "function": "functions",
        "final": "finalprocs",
        "bound": "boundprocs",
        "modproc": "modprocs",
        "common": "common",
    }

    def extend_or_append(search, name, item, value):
        """Helper function to add to search space for converting links
        """
        if name == "constructor":
            constructor = getattr(item, "constructor", False)
            if constructor:
                search.append(constructor)
        else:
            search.extend(getattr(item, value, []))

    def find_in_list(name, searchlist):
        for obj in searchlist:
            if name == obj.name.lower():
                return obj
        return None

    def get_searchlist(type_name, type_lookup, parent):
        searchlist = []
        if not type_name:
            for value in type_lookup.values():
                searchlist.extend(getattr(parent, value))
        else:
            searchlist.extend(getattr(parent, type_lookup[type_name.lower()]))
        return searchlist

    def convert_link(match):
        # [name,obj,subname,subobj]

        ERR = "Warning: Could not substitute link {}. {}"

        try:
            searchlist = get_searchlist(match.group(2), LINK_TYPES, project)
        except KeyError:
            print(
                ERR.format(
                    match.group(),
                    'Unrecognized classification "{}".'.format(match.group(2)),
                )
            )
            return match.group()

        item = find_in_list(match.group(1).lower(), searchlist)
        if item is None:
            print(ERR.format(match.group(), '"{}" not found.'.format(match.group(1))))
            return "<a>{}</a>".format(match.group(1))

        url = item.get_url()
        name = item.name

        if match.group(3):
            searchlist = []
            if not match.group(4):
                for value in SUBLINK_TYPES.values():
                    extend_or_append(searchlist, value, item, value)
            else:
                try:
                    subobject_type = SUBLINK_TYPES[match.group(4).lower()]
                except KeyError:
                    print(
                        ERR.format(
                            match.group(),
                            'Unrecognized classification "{}".'.format(match.group(4)),
                        )
                    )
                    return match.group()

                if not hasattr(item, subobject_type):
                    print(
                        ERR.format(
                            match.group(),
                            '"{}" can not be contained in "{}"'.format(
                                match.group(4), item.obj
                            ),
                        )
                    )
                    return match.group()

                extend_or_append(
                    searchlist, match.group(4).lower(), item, subobject_type,
                )

            # The other child items are anchors within the parent
            # page. But both namelists and input parameters are
            # anchors on the namelist page, and we can link to
            # either. At this point, we know we're linking to a child
            # (input parameter), so we need to recover the base URL
            url = url.split("#")[0]

            item = find_in_list(match.group(3).lower(), searchlist)
            if item is not None:
                url = url + "#" + item.anchor
                name = item.name
            else:
                print(
                    ERR.format(
                        match.group(),
                        '"{0}" not found in "{1}", linking to page for "{1}" instead.'.format(
                            match.group(3), name
                        ),
                    )
                )
        return '<a href="{}">{}</a>'.format(url, name)

    # Get information from links (need to build an RE)
    string = ford.utils.LINK_RE.sub(convert_link, string)
    return string


if __name__ == "__main__":
    # Need to monkey-patch some of our stuff with FORD-specific bits
    Namelist.get_url = get_namelist_url
    NamelistVariable.get_url = get_namelist_variable_url

    # Here we inject our patched versions into FORD
    ford.main = interceptor_main
    ford.utils.sub_links = interceptor_sub_links

    # Now we can run FORD as normal
    ford.run()
