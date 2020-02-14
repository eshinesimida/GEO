# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:19:00 2020

@author: a4546
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

TPRandFPR = pd.read_csv('AUC_spia_GSE6030.csv')

#TPRandFPR.sort_values('FPR',inplace=True)

from sklearn.metrics import auc


AUC= auc(TPRandFPR['FPR'],TPRandFPR['TPR'])

plt.scatter(x=TPRandFPR['FPR'],y=TPRandFPR['TPR'],label='(FPR,TPR)',color='k')
plt.plot(TPRandFPR['FPR'], TPRandFPR['TPR'], 'k',label='AUC = %0.2f'% AUC)
plt.legend(loc='lower right')

#plt.title('Receiver Operating Characteristic')
plt.title('ROC of SPIA_GSE6030')
plt.plot([(0,0),(1,1)],'r--')
plt.xlim([-0.01,1.01])
plt.ylim([-0.01,01.01])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()
