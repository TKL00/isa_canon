def my_function(a):
    for i in range(len(a)):
        a[0] = "pear"
        print(a[i])

my_list = ["apple", "pear"]

my_function(my_list)
print(my_list)