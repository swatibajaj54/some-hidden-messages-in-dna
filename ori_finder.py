
#with open('test1.txt') as input_file:
def get_skew_diag_data(genome):
    output = 0
    result = []
    # print(output, end=' ')
    result.append(output)
    for i in genome:
        if i == 'C':
            output = output-1
            # print(output, end=' ')
            result.append(output)
        elif i == 'G':
            output = output+1
            # print(output, end=' ')
            result.append(output)
        else:
            # print(output, end=' ')
            result.append(output)
    return result


