import csv

#get list of number of clones in each bootstrap
def num_clones(data):
    num_clones_list = []
    for row in data:
        l = len(row)
        num_clones_list.append(l)
    return num_clones_list

#make list of number of single member clones
def count_ones(data):
    ones_list = []
    for row in data:
        ones = 0
        for item in row:
            if item == 1:
                ones += 1
        ones_list.append(ones)
    return ones_list

#bin the clones by number of sequences
def buckets(data):
    buckets = [5,10,100,250,float("inf")]
    count_list = []
    for row in data:
        counts = [0,0,0,0,0]
        for item in row:
            for i in range(len(buckets)):
                if item < buckets[i]:
                    counts[i] += 1
                    break
        count_list.append(counts)
    return count_list

#get mean number of clones in each bucket
def get_mean(bucket_list):
    num_cols = len(bucket_list[0])
    mean = []
    for i in range(num_cols):
        sum_i = []
        for row in bucket_list:
            sum_i.append(row[i])
        m = sum(sum_i)/len(sum_i)
        mean.append(m)
    return mean
        
#read in data
def get_data(csv_file):
    with open(csv_file, 'rb') as f:
        reader = csv.reader(f)
        data = []
        for line in reader:
            rv = []
            for item in line:
                item=int(item)
                rv.append(item)
            data.append(rv)
    return data

if __name__ == "__main__":

    data = get_data('rare_007.csv')
    buckets_007 = buckets(data)
    mean_buckets_007 = get_mean(buckets_007)
    
    data_012 = get_data('rare_012.csv')
    buckets_012 = buckets(data_012)
    mean_buckets_012 = get_mean(buckets_012)

    data_011 = get_data('rare_011.csv')
    buckets_011 = buckets(data_011)
    mean_buckets_011 = get_mean(buckets_011)
