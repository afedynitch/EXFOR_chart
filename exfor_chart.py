from exforchart import set_debug_level


if __name__ == "__main__":

    from optparse import OptionParser

    usage = """usage: %prog [options] args"""
    parser = OptionParser(usage=usage)

    parser.add_option(
        "-i",
        "--interactive",
        action="store_true",
        dest="interact",
        help="Interative panel mode.",
    )
    parser.add_option(
        "-c",
        "--chart",
        action="store_true",
        dest="static_chart",
        help="Static chart mode for plotting.",
    )
    parser.add_option(
        "-s",
        "--save",
        action="store_true",
        dest="save_chart",
        help="Save chart as pdf and png.",
    )
    parser.add_option(
        "-l",
        "--logo",
        action="store_true",
        dest="show_logo",
        help="Show logo on chart.",
    )
    parser.add_option(
        "-d",
        "--debug",
        type=int,
        dest="debug_level",
        default=0,
        help="Level of debug printout.",
    )

    (options, args) = parser.parse_args()

    set_debug_level(options.debug_level)

    from exforchart.chart import start_interactive, start_static

    if options.interact and options.static_chart:
        raise Exception("exfor_chart(): Choose either interactive or chart mode.")
    if options.static_chart:
        start_static(options.show_logo, options.save_chart)
        exit(0)
    if options.interact:
        start_interactive(options.show_logo)
        exit(0)