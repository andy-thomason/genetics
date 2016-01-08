import sys

# simple (and horribly inefficient) forward bwt
def bwt(str):
    sorter = sorted((str[i:], i) for i in range(len(str)))
    return ''.join([str[i-1] for _, i in sorter])

# simple inverse bwt
# based on http://codereview.stackexchange.com/questions/21058/increasing-speed-of-bwt-inverse
def ibwt(str):
    sorter = sorted((char, index) for index, char in enumerate(str))
    print(sorter)
    k = sorter[0][1]
    result = []
    for i in range(len(str)):
        char, k = sorter[k]
        result.append(char)
    return ''.join(result)

bwt_str = bwt(sys.argv[1] + '$')

ibwt_str = ibwt(bwt_str)

print(bwt_str, ibwt_str)
