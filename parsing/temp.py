lst = ['1-585988', '2702782-2746290', '12954385-13000000', '13000001-13004384', '16799164-16849163', '29552234-29553835', '121757469-122000000', '122000001-122503247', '124785433-124978326', '125013061-126000000', '126000001-127000000', '127000001-128000000', '128000001-129000000', '129000001-130000000', '130000001-131000000', '131000001-132000000', '132000001-133000000', '133000001-134000000', '134000001-135000000', '135000001-136000000', '136000001-137000000', '137000001-138000000', '138000001-139000000', '139000001-140000000', '140000001-141000000', '141000001-142000000', '142000001-143000000', '143000001-143184587', '223558936-223608935', '228558365-228608364', '248946423-248956422']
cleaned = []
tmp_start = 0
tmp_end = 0
flag = False
for i in range(0,len(lst)):
	first,last = lst[i].split("-")
	if ((int(last)%1000000) == 0) and not flag:
		tmp_start = first
		flag = True
		continue
	elif flag and ((int(last)%1000000)) == 0:
		continue
	elif flag:
		cleaned.append(tmp_start+"-"+last)
		flag = False
		continue
	else:
		cleaned.append(first+"-"+last)

length = []
for i in cleaned:
	first,last = i.split("-")
	length.append(int(last)-int(first))

print(cleaned)
print(length)
idx = length.index(max(length)) #gets centromere based on longest stretch of "N"
print(cleaned[idx])

# first telomere in cleaned[0]
# centromere in cleaned[idx]
# second telomere in cleaned[-1] (last in list)

file = open("GRCH38_chr_sites.txt", 'a')
short = str(int(cleaned[0].split("-")[1])+1) + "-" + str(int(cleaned[idx].split("-")[0])-1) + "\t"
long = str(int(cleaned[idx].split("-")[1])+1) + "-" + str(int(cleaned[-1].split("-")[0])-1) + "\t"
telomere1 = cleaned[0] + "\t"
centromere = cleaned[idx] + "\t"
telomere2 = cleaned[-1] + "\n"


file.write("1\t"+"NC_000001\t"+telomere1+short+centromere+long+telomere2)



print(int(cleaned[0].split("-")[1])+1)
print(int(cleaned[idx].split("-")[0])-1)
print(int(cleaned[idx].split("-")[1])+1)
print(int(cleaned[-1].split("-")[0])-1)


chrom = "NC_000001.11"
num = chrom.split(".")[0].split("0")[-1]
print(num)
