# This code is written by Weiran Wang
# For any questions or problems, please contact the author of the code at (weiran.wang@epfl.ch)

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import csv
import os


def visualization():
    '''
    This method will read the TFA and LUDB information from .csv files and visualize the results
    The abscissa is the Number of servers
    The ordinates are the Delays of the flow of interest
    '''

    # read topology_id, TFA and LUDB from .csv files
    filepath = '/Users/wangweiran/Desktop/SemesterProject/EPFL_Network_Calculus_Semester_Project/source_sink_tandem_analysis/'
    file_number = len([name for name in os.listdir(filepath) if os.path.isfile(os.path.join(filepath, name)) and 'analysis' in name])
    print("file_number : ", file_number)
    
    X = []
    TFA = []
    LUDB = []

    # there is some problems with topology 0
    # get the topology information
    for topology_id in range(1, file_number):
        filename = 'source_sink_tandem' + str(topology_id) + '_analysis.csv'
        with open(filepath+filename, 'r') as csvAnalysis:
            reader = csvAnalysis.readlines()[-1]
            X.append(topology_id+1)
            TFA.append(float(reader.split(',')[1]))
            LUDB.append(float(reader.split(',')[2].replace('\n','')))

    # visualize the plot
    fig, ax = plt.subplots()
    ax.plot(X, TFA, marker='*', label='TFA')
    ax.plot(X, LUDB, marker='o', label='LUDB')
    plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('Number of servers')
    ax.set_ylabel('Delay of the f.o.i')
    ax.set_title('Source Sink Tandem Analysis')
    ax.legend()

    plt.show()


def main():
    visualization()


if __name__ == "__main__":
    main()
