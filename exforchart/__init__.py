import pathlib

data_dir = pathlib.Path("data").absolute()  # modified from above
debug_level = 0


def set_debug_level(level):
    global debug_level
    print("Setting debug level to", level)
    debug_level = level
