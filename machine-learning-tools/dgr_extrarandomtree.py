# work in progress: OOP

import numpy as np
import dust_gal_analyzer.DustyGalaxy as dg
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesRegressor as ETR
from sklearn.cluster import DBSCAN


methodlabel = 'ExtraRandomTreesRegression'
### load data
gal = dg(glist=np.load('gal_snapshot_'+'151'+'.npz'))
Md = gal.get('dust_mass')
Mg = gal.get('gas_mass')
Ms = gal.get('star_mass')
sfr = gal.get('SFR')
dgr = Mg/gal.get('dust_mass')
Z  = np.log10(gal.get('gas_Z')/0.0134)
r  = np.log10(gal.get('gas_radii'))
Sigma_g = Mg/(np.pi*r**2)

NN = len(Md)
KK = 6 # number of galaxy properties used to fit the DGR - galaxy properties relation
label = [r'$Z$',r'$\tau_{\rm depletion}$',r'$M_*$',r'$R$',r'$f_g$',r'$\Sigma_g$']
X = np.zeros((NN,KK+2))
X[:,0] = Z
X[:,1] = np.log10(Mg/sfr) # depletion time scale
X[:,2] = np.log10(Ms)
X[:,3] = np.log10(r)
X[:,4] = Mg/(Mg+Ms) # gas mass to baryon mass ratio
X[:,5] = np.log10(Sigma_g)
X[:,KK] = np.log10(dgr)
X[:,KK+1] = np.log10(sfr/Ms) #sSFR, used for clustering

### clean up data and prepare vectors needed for fitting
X = X[~np.isnan(X).any(axis=1)]
X = X[~np.isinf(X).any(axis=1)]
X = X[np.lexsort(np.rot90(X))] # sort the data in lexicographical order
# ND points for clustering (seperate "main sequence" from quenched galaxies)
NN = len(X[:,0])
XC = np.zeros((NN,3))
XC[:,0] = X[:,0]
XC[:,1] = X[:,KK+1]
XC[:,2] = X[:,KK]
# (X,Y) for regression
Y = X[:,KK]
X = X[:,:KK]

### clustering
clt = DBSCAN(eps=0.065,min_samples=30)
clt.fit(XC)

# clustering plot
plt.figure()
plt.plot(X[:,0],Y,'co',alpha=0.3)
for i in np.unique(clt.labels_):
	print i,len(Y[clt.labels_==i])*1.0/NN
	plt.plot(X[clt.labels_==i,0],Y[clt.labels_==i],'o',alpha=0.3,label=str(i))
plt.legend()
plt.xlabel(r'$12 + \log (O/H)$')
plt.ylabel('G/D')
plt.savefig('clustering.png',eps=300)

# only fit "main sequence"
X = X[clt.labels_==0]
Y = Y[clt.labels_==0]
NN = len(Y)


# split training and c.v. sets
ID = np.unique((np.random.random(NN*7/10)*(NN-1)).astype(int))
ID_cv = np.arange(NN)[np.where(np.array(map(lambda x: x not in ID, np.arange(NN)))==True)]
X_cv = X[ID_cv]  # ISSUE: completely randomized ID? I don't think the order of the galaxy list is random
Y_cv = Y[ID_cv]
X = X[ID]
Y = Y[ID]


def mse(X,X_pred):
	return np.mean((X-X_pred)**2)

####### regression, and make prediction
cost_train = []
cost_cv = []
depth = range(20,21)
importance = np.zeros((KK,len(depth)))
j = 0
for i in depth: # here depth = 20
	regr = ETR(n_estimators=400,criterion='mse',max_depth=i)
	regr.fit(X,Y)
	Y_pred = regr.predict(X)
	cost_train.append(np.sqrt(mse(Y,Y_pred)))
	Y_pred = regr.predict(X_cv)
	cost_cv.append(np.sqrt(mse(Y_cv,Y_pred)))
	importance[:,j] = regr.feature_importances_
	print "depth",i,"factor importance:",regr.feature_importances_
	j += 1

# plot cost
plt.figure(figsize=(10*2,8*1))

plt.subplot(1,2,1)
plt.plot(-Y_cv,-Y_pred,'co',markersize=0.5,alpha = 0.8)
plt.plot(np.linspace(-5,-1,40),np.linspace(-5,-1,40),'k--')
plt.xlabel(r'$\log\ {\rm DGR}_{\rm c.v.}$')
plt.ylabel(r'$\log\ {\rm DGR}_{\rm predict}$')

plt.subplot(1,2,2)
plt.plot(depth,cost_train,'-',label='training')
print 'sqrt(MSE) =',cost_cv[-1]
plt.plot(depth,cost_cv,'--',label='cross validation')
plt.xlabel('Depth')
plt.ylabel(r'$\sqrt{\rm MSE}\ ({\rm dex})$')
plt.legend()

plt.savefig('cv_'+methodlabel+'.png',dpi=300)
plt.close()

# plot importance
plt.figure(figsize=(10,8))
for i in range(KK):
	plt.plot(depth,importance[i],label=label[i])
plt.xlabel('Depth')
plt.ylabel('Importance (%)')
plt.legend()
plt.savefig('importance_'+methodlabel+'.png',dpi=300)
plt.close()
