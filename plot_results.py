import csv
import matplotlib.pyplot as plt

class ResultEntry:
    def __init__(self):
       self.xAxis = []
       self.stTimes = []
       self.parTimes = []
       self.stTransTimes = []
       self.parTransTimes = []
    def addPoint(self,res):
       self.xAxis.append(int(res[0]))
       self.stTimes.append(float(res[2]))
       self.parTimes.append(float(res[3]))
       self.stTransTimes.append(float(res[4]))
       self.parTransTimes.append(float(res[5]))
def split_by_cores(results):
    core_ns = set(())
    split_res = {}
    for res in results:
        core_ns.add(int(res[1]))
    for core_n in core_ns:
        res_entry = ResultEntry()
        for res in filter(lambda x: int(x[1])==core_n, results):
            res_entry.addPoint(res)
        split_res[core_n] = res_entry
    return split_res


def split_by_compilers(results):
    compiler_ids = set(())
    split_res = {}
    for res in results:
        compiler_ids.add(res[0])
    for cid in compiler_ids:
        filtered_res = []
        for res in filter(lambda x: x[0]==cid, results):
            filtered_res.append(res[1:])
        split_res[cid] = split_by_cores(filtered_res)
    return split_res

def plot_res(results):
    index = 1
    for cid, cres in results.items():
        for ncores, res in cres.items():
            plt.subplot(len(results), len(cres), index)
            index += 1
            plt.plot(res.xAxis, res.stTimes, label='st')
            plt.plot(res.xAxis, res.parTimes, label='par')
            plt.plot(res.xAxis, res.stTransTimes, label='st trans')
            plt.plot(res.xAxis, res.parTransTimes, label='par trans')
            plt.xlabel('matrix size')
            plt.ylabel('time, ns')
            plt.title(cid + ': ' + str(ncores) + ' cores')
            plt.grid()
            leg = plt.legend(loc='upper center')
    plt.show()
with open('result.csv') as resultfile:
    result_data = csv.reader(resultfile, delimiter=';')
    result_list = []
    for row in result_data:
        result_list.append(row)
    per_compiler_results = split_by_compilers(result_list)
    plot_res(per_compiler_results)
