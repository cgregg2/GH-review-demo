# Read this short program and try to predict what it does.
# Run it: how accurate was your prediction?
# Refactor the program to make it more readable. Remember to run it after each change to ensure its behavior hasn’t changed.
# Compare your rewrite with your neighbor’s. What did you do the same? What did you do differently, and why?

n = 10
s = 'et cetera'
print(s)
i = 0
while i < n:
    # print('at', j)
    new = ''
    for j in range(len(s)):
        left = j-1
        right = (j+1)%len(s)
        if s[left]==s[right]: new += '-'
        else: new += '*'
    s=''.join(new)
    print(s)
    i += 1

