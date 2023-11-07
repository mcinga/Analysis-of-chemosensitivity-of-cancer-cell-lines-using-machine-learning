# -*- coding: utf-8 -*-


from sklearn.linear_model import ElasticNet
from sklearn.neighbors import KNeighborsRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import pandas as pd
import numpy as np

#==============================================================================
df01=pd.read_csv("GDSC_01_CRISPR_reg5000.csv")
# Scaling and centering the data
df01.loc[:, df01.dtypes == 'float64'] = df01.loc[:, df01.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())
# Changing categorical variables into numericals
df01 = pd.get_dummies(df01)


X = df01.drop("DOSE_RESPONSE", axis=1)
y = df01["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df01)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg01.csv", index=False)

#==============================================================================

df02=pd.read_csv("GDSC_02_CRISPR_reg5000_2.csv")

df02.loc[:, df02.dtypes == 'float64'] = df02.loc[:, df02.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df02 = pd.get_dummies(df02)


X = df02.drop("DOSE_RESPONSE", axis=1)
y = df02["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df02)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg02.csv", index=False)

#==============================================================================

df03=pd.read_csv("GDSC_03_CRISPR_reg5000.csv")

df03.loc[:, df03.dtypes == 'float64'] = df03.loc[:, df03.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df03 = pd.get_dummies(df02)


X = df03.drop("DOSE_RESPONSE", axis=1)
y = df03["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df03)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg03.csv", index=False)

#=============================================================================
df04=pd.read_csv("GDSC_04_CRISPR_reg5000.csv")

df04.loc[:, df04.dtypes == 'float64'] = df04.loc[:, df04.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df04 = pd.get_dummies(df04)


X = df04.drop("DOSE_RESPONSE", axis=1)
y = df04["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df04)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg04.csv", index=False)

#==============================================================================
df05=pd.read_csv("GDSC_05_CRISPR_reg5000.csv")

df05.loc[:, df05.dtypes == 'float64'] = df05.loc[:, df05.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df05 = pd.get_dummies(df05)


X = df05.drop("DOSE_RESPONSE", axis=1)
y = df05["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df05)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg05.csv", index=False)

#=============================================================================
df06=pd.read_csv("GDSC_06_CRISPR_reg5000.csv")

df06.loc[:, df06.dtypes == 'float64'] = df06.loc[:, df06.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df06 = pd.get_dummies(df06)


X = df06.drop("DOSE_RESPONSE", axis=1)
y = df06["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df06)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg06.csv", index=False)

#==============================================================================
df07=pd.read_csv("GDSC_07_CRISPR_reg5000.csv")

df07.loc[:, df07.dtypes == 'float64'] = df07.loc[:, df07.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df07 = pd.get_dummies(df07)


X = df07.drop("DOSE_RESPONSE", axis=1)
y = df07["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df07)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg07.csv", index=False)

#==============================================================================
df08=pd.read_csv("GDSC_08_CRISPR_reg5000.csv")

df08.loc[:, df08.dtypes == 'float64'] = df08.loc[:, df08.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df08 = pd.get_dummies(df08)


X = df08.drop("DOSE_RESPONSE", axis=1)
y = df08["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df08)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg08.csv", index=False)

#==============================================================================
df09=pd.read_csv("GDSC_09_CRISPR_reg5000.csv")

df09.loc[:, df09.dtypes == 'float64'] = df09.loc[:, df02.dtypes == 'float64'].apply(lambda x: (x - x.mean()) / x.std())

df09 = pd.get_dummies(df02)


X = df09.drop("DOSE_RESPONSE", axis=1)
y = df09["DOSE_RESPONSE"]
inTraining = np.random.rand(len(df09)) < 0.7
training = X[inTraining]
testing = X[~inTraining]
y_train = y[inTraining]
y_test = y[~inTraining]


models = {"ElasticNet": ElasticNet(alpha=1.0, l1_ratio=0.5),
          "KNN": KNeighborsRegressor(n_neighbors=5),
          "GradientBoosting": GradientBoostingRegressor(n_estimators=100, max_depth=3),
          "RandomForest": RandomForestRegressor(n_estimators=100, max_depth=10)}


results_table = pd.DataFrame(columns=["models", "RMSE", "R2", "MAE"])

#Loop through the models
for model_name, model in models.items():
  #train the models
    model.fit(training, y_train)
   #make predictions
    y_pred = model.predict(testing)
   #evaluate the performane of the model
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    #store the results into the dataframe 
    results_table = results_table.append({"models": model_name, "RMSE": rmse, "R2": r2, "MAE": mae}, ignore_index=True)

#store the results into a csv file to open in excel.
results_table.to_csv("results_table_reg09.csv", index=False)

