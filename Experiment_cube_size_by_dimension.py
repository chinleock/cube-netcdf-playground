'''
This script is to test 
    1) with the same size of data, but specify different size of dimension, 
       how does the nc file size being changed according to dime size.
    2) also this script is to test the unlimited dim is appendable.

Output Files: 
    File 1: testTime_largerSize.nc   : 10 times larger in time dimension
    File 2: testTime_exactSize.nc    : Exact size in time dimension
    File 3: testTime_unlimitedSize.nc: Set the time dimension as unlimited 

Test 1) 
    Expected result in different time dimension size should be:
        testTime_largerSize.ncFile    size: 194.433299MB
        testTime_exactSize.ncFile     size: 19.451699MB
        testTime_unlimitedSize.ncFile size: 15.342867MB
    
    Set the time dimension to unlimited can save more space

Test 2)
    

Created by Leo Chin, 1/16/2015
'''

################################### Import ####################################
import netCDF4 as NET
import numpy as NUM
from os import path as PATH

################################### Global ####################################
ncFileDir = r'C:\Temp'
testFile1 = PATH.join(ncFileDir, "testTime_largerSize.nc")
testFile2 = PATH.join(ncFileDir, "testTime_exactSize.nc")
testFile3 = PATH.join(ncFileDir, "testTime_unlimitedSize.nc")

timeData = NUM.arange(0, 300)       # numpy array size of 300
xData = NUM.arange(10, 100)         # numpy array size of 90
yData = NUM.arange(100, 10, -1)     # numpy array size of 90
countData = NUM.arange(0, 300 * 90 * 90).reshape(300, 90, 90)
# numpy array, 3-Dim, shape of 300 by 90 by 90

addTimeData = NUM.arange(300, 400)  # numpy array size of 100
addCountData = NUM.arange(0, 100 * 90 * 90).reshape(100, 90, 90)
# numpy array, 3-Dim, shape of 100 by 90 by 90

################################# Functions ###################################

#### Create Dimension and Its Variables ####
def createDim(dataset, timeShape = None):
    dataset.createDimension('time', timeShape)
    dataset.createDimension('x', 90)
    dataset.createDimension('y', 90)
    
    time = dataset.createVariable('time', 'f8', ('time'))
    x = dataset.createVariable('x', 'f8', ('x'))
    y = dataset.createVariable('y', 'f8', ('y'))
    
    time[:300] = timeData
    x[:] = xData
    y[:] = yData

#### Create Count Variable ####
def createVariable(dataset):
    var = dataset.createVariable('COUNT', 'f8', ('time', 'y', 'x'))
    var[:300, :, :] = countData

#### Append Time Variable and Also Count Variable ####
def appendTime(dataset):
    #### Change Time Variable Size to Extend 100 More ####
    newTimeData = NUM.concatenate((timeData, addTimeData))
    dataset.variables['time'][:] = newTimeData
    
    #### Change Count Variable Expand in Time Dimension to 100 More ####
    newCountData = NUM.concatenate((countData, addCountData), axis = 0)
    dataset.variables['COUNT'][:] = newCountData

####  Supplementary Function: Get File Size Return MB ####
def getFileSize(filePath):
    return str(PATH.getsize(filePath) / 1.0e6) + "MB"

###################################### Main ######################################
if __name__ == '__main__':
    #### Create testFile1 ####
    dataset1 = NET.Dataset(testFile1, 'w')
    createDim(dataset1, 3000)
    createVariable(dataset1)
    size = getFileSize(testFile1)
    print("Done " + testFile1 + "File size: " + size)
    dataset1.close()
    
    #### Create testFile2 ####
    dataset2 = NET.Dataset(testFile2, 'w')
    createDim(dataset2, 300)
    createVariable(dataset2)
    size = getFileSize(testFile2)
    print("Done " + testFile2 + "File size: " + size)
    dataset2.close()

    #### Create testFile3 ####
    dataset3 = NET.Dataset(testFile3, 'w')
    createDim(dataset3)
    createVariable(dataset3)
    size = getFileSize(testFile3)
    print("Done " + testFile3 + "File size: " + size)
    appendTime(dataset3)
    size = getFileSize(testFile3)
    print("Done appending " + testFile3 + "File size: " + size)
    dataset3.close()
