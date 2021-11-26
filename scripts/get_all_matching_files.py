#!/usr/bin/env python

def get_all_matching_files(paths = None, patterns = None, absolute = False, exclude=None):
    if paths is None:
        paths = ['./']
    if patterns is None:
        patterns = ['*']
    if exclude is None:
        exclude = []
        
    if not isinstance(paths, list):
        paths = [paths]
    if not isinstance(patterns, list):
        patterns = [patterns]
    if not isinstance(exclude, list):
        exclude = [exclude]

    # Use recursion when number of paths/patterns
    # is bigger than one
    if len(patterns) <= 0:
        return []
    elif len(patterns) > 1:
        result = get_all_matching_files(paths, patterns[0], absolute, exclude)
        result.extend(get_all_matching_files(paths, patterns[1:], absolute, exclude))
        return result
    else:
        if len(paths) <= 0:
            return []
        elif len(paths) > 1:
            result = get_all_matching_files(paths[0], patterns[0], absolute, exclude)
            result.extend(get_all_matching_files(paths[1:], patterns[0], absolute, exclude))
            return result

    # Here we deal with the case when we have exactly one path
    # and pattern to match
    
    import os
    import fnmatch

    result = []
    for root, dirname, filenames in os.walk(paths[0], topdown=True):
        dirname[:] = [d for d in dirname if os.path.join(root, d) not in exclude]
        for filename in fnmatch.filter(filenames, patterns[0]):
            result.append(os.path.join(root, filename))

    if absolute:
        result = [os.path.abspath(x) for x in result]
        
    return result

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
        description='Recursively find files matching various patterns. Acts a bit like a simple but portable alternative to find.')

    parser.add_argument("--path", action = "append", help =
                        "Specify a path to search recursively.")
    
    parser.add_argument("--pattern", action = "append", help =
                        "Specify a pattern to search recursively.")

    parser.add_argument("--absolute", action = "store_true", help =
                        "If set then report absolute paths.")
    
    parser.add_argument("--exclude", action="append", help="Paths to exclude from search")

    args = parser.parse_args()

    result = get_all_matching_files(paths=args.path, patterns = args.pattern, absolute = args.absolute, exclude=args.exclude)

    for f in result:
        print(f)
