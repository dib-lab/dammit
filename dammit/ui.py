def add_header(msg, level):
    s = ''
    if level == 0:
        symbol = '='
        return '{0}\n{1}\n{2}\n'.format(symbol * 40,
                                      msg,
                                      symbol * 40)
    elif level == 1:
        symbol = '~'
        return '{0}\n{1}\n{2}'.format(symbol * 30,
                                      msg,
                                      symbol * 30)
    else:
        symbol = '---'
        return '\n{0} {1}\n'.format(symbol,
                                    msg)


def print_header(msg, level):
    '''Standardize output headers for submodules.

    This doesn't need to be logged, but it's nice for
    the user.
    '''
    print(add_header(msg, level), file=sys.stderr)


