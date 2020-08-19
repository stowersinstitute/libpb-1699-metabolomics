from numpy import clip, sign, column_stack, outer
from numpy.random import permutation, seed
from sklearn.metrics import roc_curve, roc_auc_score
from pyopls import OPLS
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import cross_val_predict, LeaveOneOut
from sklearn.metrics import r2_score, accuracy_score
from pandas import DataFrame
from seaborn import heatmap, scatterplot
from scipy.stats import ttest_ind_from_stats, norm
#from scipy.stats import norm
from math import nan

def compute_metrics(pls, data, categories, valcategories):
    y_pred = cross_val_predict(pls, data, categories, cv=LeaveOneOut())
    q_squared = r2_score(valcategories, y_pred)
    dq_squared = r2_score(valcategories, clip(y_pred, -1, 1))
    accuracy = accuracy_score(valcategories, sign(y_pred))

    return y_pred, q_squared, dq_squared, accuracy

def run_opls(data, categories, category_names, population, tissue, n_components, n_p_val_iter):
    seed(1234567)
    opls = OPLS(n_components)
    Z = opls.fit_transform(data, categories)
    # to get original data back:
    #Z += outer(opls.T_ortho_[:,0], opls.P_ortho_[:,0])

    pls = PLSRegression(1)
    y_pred, q_squared, dq_squared, accuracy = compute_metrics(pls, data, categories, categories)

    processed_y_pred, processed_q_squared, processed_dq_squared, processed_accuracy = compute_metrics(pls, Z, categories, categories)

    #r2_X = opls.score(data)

    pls.fit(Z, categories)

    # compute p-value
    # TODO: validate
    dq2s = []
    seed(1234567)
    if n_p_val_iter > 0:
        for k in range(n_p_val_iter):
            pcategories = permutation(categories)
            o = OPLS(n_components)
            Z = opls.fit_transform(data, pcategories)
            p = PLSRegression(1)
            _, _, dq2, _ = compute_metrics(p, Z, pcategories, categories)
            dq2s.append(dq2)
        seed(1234567)
        dq2s = DataFrame(data=dq2s)
        #print(dq2s)
        dq2_mean = float(dq2s.mean())
        dq2_std = float(dq2s.std())
        p = 2.*norm.sf(processed_dq_squared, dq2_mean, dq2_std)
        if p > 1.:
            p = 2.*norm.sf(dq2_mean+(dq2_mean-processed_dq_squared), dq2_mean, dq2_std)
            assert p <= 1.
        #print(dq2_mean, dq2_std)
    else:
        p = nan

    #t,p = ttest_ind_from_stats(dq2_mean, dq2_std, len(dq2s.index), processed_dq_squared, dq2s.std(), len(dq2s.index))
    #print(f'{float(p):.2e}')

    return (pls, opls, Z, processed_q_squared, processed_dq_squared, float(p), processed_accuracy, y_pred, processed_y_pred)

def plot_opls(pls, opls, data, categories, category_names, population, tissue, n_components, y_pred_pre, y_pred, q2, dq2, p, acc, scores_ax=None, loadings_ax=None, roc_ax=None):
    fpr, tpr, thresholds = roc_curve(categories, y_pred_pre)
    roc_auc = roc_auc_score(categories, y_pred_pre)
    proc_fpr, proc_tpr, proc_thresholds = roc_curve(categories, y_pred)
    proc_roc_auc = roc_auc_score(categories, y_pred)

    df = DataFrame(column_stack([pls.x_scores_, opls.T_ortho_[:, 0]]),
          index=data.index, columns=['t', 't_ortho'])
    # plot scores
    if scores_ax is not None:
        for u,c in zip([1,-1], category_names):
            scores_ax.scatter(df[categories==u]['t'], df[categories==u]['t_ortho'], label=c)
        scores_ax.set_xlabel('PLS score (t)')
        scores_ax.set_ylabel('PLS score (t orthogonal)')
        scores_ax.set_title(f'Scores for {population}, {tissue}\n$Q^2={q2:.2f}$, $DQ^2={dq2:.2f}$, $p={p:.2e}$, acc.={acc:.2f}')
        scores_ax.legend(loc='best')

    # plot loadings
    if loadings_ax is not None:
        loadings = DataFrame(pls.x_loadings_, index=data.columns.rename(''))
        scatterplot(data=loadings, ax=loadings_ax, legend=False)
        loadings_ax.set_title(f'Loading for {population}, {tissue}')
        #https://stackoverflow.com/questions/51898101/how-do-i-stagger-or-offset-x-axis-labels-in-matplotlib
        for tick in loadings_ax.xaxis.get_major_ticks()[1::2]:
            tick.set_pad(15)

    if roc_ax is not None:
        roc_ax.plot(fpr, tpr, lw=2, label=f'PLS (AUC={roc_auc:.4f})')
        roc_ax.plot(proc_fpr, proc_tpr, lw=2, label=f'OPLS (AUC={proc_roc_auc:.4f})')
        roc_ax.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
        roc_ax.set_xlabel('False Positive Rate')
        roc_ax.set_ylabel('True Positive Rate')
        roc_ax.set_title('ROC Curve')
        roc_ax.legend(loc='lower right')
