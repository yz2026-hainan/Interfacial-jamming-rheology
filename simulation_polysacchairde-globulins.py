#!usrbinenv python
# -- coding utf-8 --

Rheological Analysis of Soy Protein (7S & 11S) at Different Concentrations
Individual CSV file per protein per concentration per test type.
Includes interfacial jamming analysis with raw data export and PCA coordinates.


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
import networkx as nx
from mpl_toolkits.mplot3d import Axes3D
import warnings
import os
warnings.filterwarnings('ignore')

# ============================================================================
# Nature Communications Style
# ============================================================================
plt.style.use('default')
plt.rcParams.update({
    'font.size' 11,
    'font.family' 'Arial',
    'axes.labelsize' 12,
    'axes.titlesize' 14,
    'xtick.labelsize' 10,
    'ytick.labelsize' 10,
    'legend.fontsize' 10,
    'figure.dpi' 300,
    'savefig.dpi' 300,
    'savefig.bbox' 'tight',
    'axes.linewidth' 1.2,
    'grid.linewidth' 0.5,
    'lines.linewidth' 2,
})

PROTEIN_COLORS = {'7S' '#1f77b4', '11S' '#d62728'}  # blue, red

print(=  80)
print(RHEOLOGICAL ANALYSIS OF SOY PROTEIN (7S & 11S) AT DIFFERENT CONCENTRATIONS)
print(Includes interfacial jamming analysis with full PCA data export)
print(=  80)

# ============================================================================
# Helper read CSV with multiple encoding attempts
# ============================================================================
def read_csv_with_encoding(filepath)
    Try multiple encodings to read CSV file, without chardet.
    encodings = ['utf-8', 'gbk', 'gb2312', 'gb18030', 'latin1', 'cp1252']
    for enc in encodings
        try
            df = pd.read_csv(filepath, encoding=enc)
            return df
        except
            continue
    try
        df = pd.read_csv(filepath, encoding='utf-8', errors='replace')
        return df
    except
        raise ValueError(fCould not read {filepath} with any common encoding)

# ============================================================================
# Data Loader (individual files per protein+concentration+test)
# ============================================================================
class ProteinRheologyLoader
    def __init__(self, data_path='FDatarheology_data')
        self.data_path = data_path
        self.proteins = ['7S', '11S']
        self.concentrations = [1.0, 1.5, 2.0, 7.5, 10.0]
        self.test_types = ['laos', 'frequency', 'creep', 'shear']
        self.sample_names = []
        self.data_dict = {}
        self.sample_info = pd.DataFrame()

    def load_all_data(self)
        print(n1. LOADING INDIVIDUAL FILES (Protein_Concentration_Test.csv))
        print(-  40)
        loaded_count = 0
        for protein in self.proteins
            for conc in self.concentrations
                sample_name = f{protein}_{conc}%
                self.sample_names.append(sample_name)
                self.data_dict[sample_name] = {}
                for test in self.test_types
                    filename = f{protein}_{conc}_{test}.csv
                    filepath = os.path.join(self.data_path, filename)
                    if os.path.exists(filepath)
                        try
                            df = read_csv_with_encoding(filepath)
                            self.data_dict[sample_name][test] = df
                            loaded_count += 1
                            print(f  ✓ {filename})
                        except Exception as e
                            print(f  ✗ Error reading {filename} {e})
                            self.data_dict[sample_name][test] = None
                    else
                        alt_filename = f{protein}_{conc}_{test}.csv
                        alt_path = os.path.join(self.data_path, alt_filename)
                        if os.path.exists(alt_path)
                            try
                                df = read_csv_with_encoding(alt_path)
                                self.data_dict[sample_name][test] = df
                                loaded_count += 1
                                print(f  ✓ {alt_filename})
                            except
                                print(f  ✗ Missing {filename})
                                self.data_dict[sample_name][test] = None
                        else
                            print(f  ✗ Missing {filename})
                            self.data_dict[sample_name][test] = None
        info_list = []
        for sid in self.sample_names
            parts = sid.split('_')
            protein = parts[0]
            conc_str = parts[1].replace('%', '')
            conc = float(conc_str)
            info_list.append({'sample_id' sid, 'Protein' protein, 'Concentration_%' conc})
        self.sample_info = pd.DataFrame(info_list)
        print(fn  Loaded {loaded_count} files, {len(self.sample_names)} samples.)
        return self.data_dict

# ============================================================================
# Parameter Extractor (with unit conversions)
# ============================================================================
class RheologyParameterExtractor
    def __init__(self, data_dict, sample_info)
        self.data_dict = data_dict
        self.sample_info = sample_info
        self.parameters = pd.DataFrame()

    def extract_all_parameters(self)
        print(n2. EXTRACTING RHEOLOGICAL PARAMETERS)
        print(-  40)
        param_list = []
        for sample_id, tests in self.data_dict.items()
            row_info = self.sample_info[self.sample_info['sample_id'] == sample_id]
            protein = row_info.iloc[0]['Protein'] if len(row_info)  0 else 'Unknown'
            conc = row_info.iloc[0]['Concentration_%'] if len(row_info)  0 else np.nan
            params = {'Sample_ID' sample_id, 'Protein' protein, 'Concentration_%' conc}
            if 'laos' in tests and tests['laos'] is not None
                params.update(self._extract_laos(tests['laos']))
            if 'frequency' in tests and tests['frequency'] is not None
                params.update(self._extract_frequency(tests['frequency']))
            if 'shear' in tests and tests['shear'] is not None
                params.update(self._extract_shear(tests['shear']))
            if 'creep' in tests and tests['creep'] is not None
                params.update(self._extract_creep(tests['creep']))
            param_list.append(params)
        self.parameters = pd.DataFrame(param_list)
        self._calculate_rpi()
        print(f  Extracted {len(self.parameters.columns)} parameters for {len(self.parameters)} samples)
        return self.parameters

    def _extract_laos(self, df)
        params = {}
        if 'storage_modulus(Pa)' in df.columns
            Gp = df['storage_modulus(Pa)'].values
            strain = df['strain_amplitude(%)'].values if 'strain_amplitude(%)' in df.columns else np.arange(len(Gp))
            linear_G = np.mean(Gp[min(3, len(Gp))])
            params['G_prime_linear_Pa'] = linear_G
            min_G = np.min(Gp)
            params['Payne_effect'] = (linear_G - min_G)  linear_G if linear_G  0 else 0
            target = 0.95  linear_G
            idx = np.where(Gp  target)[0]
            params['Linear_limit_strain_%'] = strain[idx[0]] if len(idx)  0 else strain[-1]
            max_G = np.max(Gp)
            params['Strain_stiffening_ratio'] = max_G  linear_G if linear_G  0 else 1
        return params

    def _extract_frequency(self, df)
        params = {}
        if 'angular_frequency(Hz)' in df.columns and 'storage_modulus(Pa)' in df.columns
            omega_hz = df['angular_frequency(Hz)'].values
            omega_rad = omega_hz  2  np.pi
            Gp = df['storage_modulus(Pa)'].values
            valid = (omega_rad  0) & (Gp  0)
            if np.sum(valid)  2
                log_omega = np.log10(omega_rad[valid])
                log_Gp = np.log10(Gp[valid])
                coeff = np.polyfit(log_omega, log_Gp, 1)
                params['G_prime_power_law_exponent'] = coeff[0]
            if 'loss_modulus(Pa)' in df.columns
                Gdp = df['loss_modulus(Pa)'].values
                tan_delta = Gdp  Gp
                params['tan_delta_avg'] = np.mean(tan_delta[np.isfinite(tan_delta)])
                for i in range(len(Gp)-1)
                    if (Gp[i] = Gdp[i]) and (Gp[i+1] = Gdp[i+1])
                        freq_i = omega_rad[i]
                        freq_ip1 = omega_rad[i+1]
                        cross = freq_i + (freq_ip1 - freq_i)  (Gp[i]-Gdp[i])  ((Gp[i]-Gdp[i]) - (Gp[i+1]-Gdp[i+1]))
                        params['Crossover_frequency_rad_s'] = cross
                        break
        return params

    def _extract_shear(self, df)
        params = {}
        if 'shear_rate(1s)' in df.columns and 'stress(mPa)' in df.columns
            gamma = df['shear_rate(1s)'].values
            sigma_mPa = df['stress(mPa)'].values
            sigma = sigma_mPa  1000.0
            params['Yield_stress_Pa'] = sigma[0] if len(sigma)  0 else np.nan
            if len(gamma)  0
                if 'viscosity(mPa.s)' in df.columns
                    visc_mPa = df['viscosity(mPa.s)'].values
                    visc = visc_mPa  1000.0
                else
                    visc = sigma  gamma
                if 10 in gamma
                    idx = np.where(gamma == 10)[0]
                    params['Viscosity_10s-1_Pas'] = visc[idx[0]] if len(idx)  0 else np.nan
                else
                    idx_nearest = np.argmin(np.abs(gamma - 10))
                    params['Viscosity_10s-1_Pas'] = visc[idx_nearest]
            try
                def hb(g, s0, K, n)
                    return s0 + K  gn
                popt, _ = curve_fit(hb, gamma, sigma, p0=[sigma[0], 1, 0.5], maxfev=5000)
                params['HB_sigma_y_Pa'] = popt[0]
                params['HB_K_Pas^n'] = popt[1]
                params['HB_n'] = popt[2]
            except
                pass
        return params

    def _extract_creep(self, df)
        params = {}
        if 'strain(%)' in df.columns and 'time(s)' in df.columns
            strain = df['strain(%)'].values
            time = df['time(s)'].values
            params['Max_creep_strain_%'] = np.max(strain)
            idx_max = np.argmax(strain)
            if idx_max  len(strain)-1
                final_strain = strain[-1]
                max_strain = strain[idx_max]
                params['Recovery_percentage'] = (1 - final_strainmax_strain)  100 if max_strain  0 else 0
            else
                params['Recovery_percentage'] = np.nan
        return params

    def _calculate_rpi(self)
        key_params = ['G_prime_linear_Pa', 'Yield_stress_Pa', 'Recovery_percentage', 'Viscosity_10s-1_Pas']
        weights = [0.3, 0.3, 0.2, 0.2]
        available = [p for p in key_params if p in self.parameters.columns]
        w_avail = [weights[i] for i, p in enumerate(key_params) if p in self.parameters.columns]
        if available
            norm_df = self.parameters[available].apply(lambda x (x - x.min())  (x.max() - x.min()))
            rpi = np.zeros(len(self.parameters))
            for col, w in zip(available, w_avail)
                rpi += norm_df[col].fillna(0)  w
            if rpi.max()  0
                rpi = rpi  rpi.max()
            self.parameters['RPI'] = rpi
        else
            self.parameters['RPI'] = 0.5

# ============================================================================
# 2D Network Mapper
# ============================================================================
class Network2DMapper
    def __init__(self, parameters_df)
        self.df = parameters_df

    def create_2d_maps(self)
        print(n3. CREATING 2D NETWORK MAPS)
        print(-  40)
        numeric_cols = self.df.select_dtypes(include=[np.number]).columns.tolist()
        exclude = ['Concentration_%', 'RPI']
        numeric_cols = [c for c in numeric_cols if c not in exclude]
        X = self.df[numeric_cols].fillna(self.df[numeric_cols].median()).values
        X_scaled = StandardScaler().fit_transform(X)
        pca = PCA(n_components=2)
        pca_result = pca.fit_transform(X_scaled)
        self.df['PC1'] = pca_result[, 0]
        self.df['PC2'] = pca_result[, 1]
        perplexity = min(5, len(X)-1)
        tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
        tsne_result = tsne.fit_transform(X_scaled)
        self.df['tSNE1'] = tsne_result[, 0]
        self.df['tSNE2'] = tsne_result[, 1]
        self._perform_clustering(X_scaled)
        self._plot_2d_maps(pca)
        return self.df

    def _perform_clustering(self, X_scaled)
        n = len(X_scaled)
        if n  3
            self.df['Cluster'] = 0
            return
        max_k = min(5, n-1)
        sil_scores = []
        for k in range(2, max_k+1)
            km = KMeans(n_clusters=k, random_state=42)
            labels = km.fit_predict(X_scaled)
            if len(np.unique(labels))  1
                sil = silhouette_score(X_scaled, labels)
                sil_scores.append(sil)
            else
                sil_scores.append(-1)
        best_k = range(2, max_k+1)[np.argmax(sil_scores)]
        km = KMeans(n_clusters=best_k, random_state=42)
        self.df['Cluster'] = km.fit_predict(X_scaled)

    def _plot_2d_maps(self, pca)
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        ax = axes[0,0]
        for prot in self.df['Protein'].unique()
            sub = self.df[self.df['Protein'] == prot]
            ax.scatter(sub['PC1'], sub['PC2'], s=150, label=prot,
                      color=PROTEIN_COLORS.get(prot, 'gray'), alpha=0.8, edgecolors='black')
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]100.1f}%)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]100.1f}%)')
        ax.set_title('a) PCA by Protein Type')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax = axes[0,1]
        sc = ax.scatter(self.df['PC1'], self.df['PC2'], c=self.df['Concentration_%'],
                       s=150, cmap='viridis', alpha=0.8, edgecolors='black')
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label('Concentration (%)')
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]100.1f}%)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]100.1f}%)')
        ax.set_title('b) PCA by Concentration')
        ax.grid(True, alpha=0.3)
        ax = axes[1,0]
        for prot in self.df['Protein'].unique()
            sub = self.df[self.df['Protein'] == prot]
            ax.scatter(sub['tSNE1'], sub['tSNE2'], s=150, label=prot,
                      color=PROTEIN_COLORS.get(prot, 'gray'), alpha=0.8, edgecolors='black')
        ax.set_xlabel('t-SNE dimension 1')
        ax.set_ylabel('t-SNE dimension 2')
        ax.set_title('c) t-SNE by Protein Type')
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax = axes[1,1]
        sc = ax.scatter(self.df['PC1'], self.df['PC2'], c=self.df['Cluster'],
                       s=150, cmap='Set2', alpha=0.8, edgecolors='black')
        for _, row in self.df.iterrows()
            ax.annotate(row['Sample_ID'], (row['PC1'], row['PC2']),
                       xytext=(5,5), textcoords='offset points', fontsize=8)
        cb = plt.colorbar(sc, ax=ax)
        cb.set_label('Cluster ID')
        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]100.1f}%)')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]100.1f}%)')
        ax.set_title('d) Clustering Results')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('2D_Network_Maps_Protein.png', dpi=300)
        print(  ✓ Saved 2D_Network_Maps_Protein.png)

# ============================================================================
# 3D Visualizations
# ============================================================================
class Network3DVisualizer
    def __init__(self, parameters_df)
        self.df = parameters_df

    def create_3d_visualizations(self)
        print(n4. CREATING 3D VISUALIZATIONS)
        print(-  40)
        fig = plt.figure(figsize=(18, 6))
        ax1 = fig.add_subplot(131, projection='3d')
        self._plot_3d_pca(ax1)
        ax2 = fig.add_subplot(132, projection='3d')
        self._plot_3d_parameter_space(ax2)
        ax3 = fig.add_subplot(133, projection='3d')
        self._plot_3d_network(ax3)
        plt.tight_layout()
        plt.savefig('3D_Visualizations_Protein.png', dpi=300)
        print(  ✓ Saved 3D_Visualizations_Protein.png)

    def _plot_3d_pca(self, ax)
        numeric_cols = self.df.select_dtypes(include=[np.number]).columns.tolist()
        exclude = ['Concentration_%', 'RPI', 'Cluster', 'PC1', 'PC2', 'tSNE1', 'tSNE2']
        numeric_cols = [c for c in numeric_cols if c not in exclude]
        X = self.df[numeric_cols].fillna(self.df[numeric_cols].median()).values
        X_scaled = StandardScaler().fit_transform(X)
        pca3 = PCA(n_components=3)
        pca3_result = pca3.fit_transform(X_scaled)
        for i, row in self.df.iterrows()
            color = PROTEIN_COLORS.get(row['Protein'], 'gray')
            ax.scatter(pca3_result[i,0], pca3_result[i,1], pca3_result[i,2],
                      color=color, s=150, alpha=0.8, edgecolors='black')
            ax.text(pca3_result[i,0], pca3_result[i,1], pca3_result[i,2],
                   row['Sample_ID'], fontsize=8)
        ax.set_xlabel(f'PC1 ({pca3.explained_variance_ratio_[0]100.1f}%)')
        ax.set_ylabel(f'PC2 ({pca3.explained_variance_ratio_[1]100.1f}%)')
        ax.set_zlabel(f'PC3 ({pca3.explained_variance_ratio_[2]100.1f}%)')
        ax.set_title('3D PCA of Protein Samples')
        ax.grid(True)

    def _plot_3d_parameter_space(self, ax)
        candidates = ['G_prime_linear_Pa', 'Yield_stress_Pa', 'Viscosity_10s-1_Pas']
        avail = [p for p in candidates if p in self.df.columns]
        if len(avail) = 3
            x = self.df[avail[0]].fillna(0)
            y = self.df[avail[1]].fillna(0)
            z = self.df[avail[2]].fillna(0)
            for i, row in self.df.iterrows()
                color = PROTEIN_COLORS.get(row['Protein'], 'gray')
                ax.scatter(x.iloc[i], y.iloc[i], z.iloc[i], color=color, s=150, alpha=0.8)
                ax.text(x.iloc[i], y.iloc[i], z.iloc[i], row['Sample_ID'], fontsize=8)
            ax.set_xlabel(avail[0])
            ax.set_ylabel(avail[1])
            ax.set_zlabel(avail[2])
            ax.set_title('3D Parameter Space')
        else
            ax.text(0.5,0.5,0.5, Insufficient parameters, ha='center')

    def _plot_3d_network(self, ax)
        numeric_cols = self.df.select_dtypes(include=[np.number]).columns.tolist()
        exclude = ['Concentration_%', 'RPI', 'Cluster', 'PC1', 'PC2', 'tSNE1', 'tSNE2']
        numeric_cols = [c for c in numeric_cols if c not in exclude]
        X = self.df[numeric_cols].fillna(self.df[numeric_cols].median()).values
        X_scaled = StandardScaler().fit_transform(X)
        from sklearn.metrics.pairwise import cosine_similarity
        sim = cosine_similarity(X_scaled)
        threshold = np.percentile(sim.flatten(), 70)
        G = nx.Graph()
        for i, row in self.df.iterrows()
            G.add_node(row['Sample_ID'], protein=row['Protein'])
        for i in range(len(self.df))
            for j in range(i+1, len(self.df))
                if sim[i,j]  threshold
                    G.add_edge(self.df.iloc[i]['Sample_ID'], self.df.iloc[j]['Sample_ID'], weight=sim[i,j])
        if len(G.nodes) == 0
            ax.text(0.5,0.5,0.5, No network edges, ha='center')
            return
        pos = nx.spring_layout(G, dim=3, seed=42)
        for edge in G.edges()
            x = [pos[edge[0]][0], pos[edge[1]][0]]
            y = [pos[edge[0]][1], pos[edge[1]][1]]
            z = [pos[edge[0]][2], pos[edge[1]][2]]
            ax.plot(x, y, z, color='gray', alpha=0.5, linewidth=1)
        for node in G.nodes()
            protein = G.nodes[node]['protein']
            color = PROTEIN_COLORS.get(protein, 'gray')
            ax.scatter(pos[node][0], pos[node][1], pos[node][2], color=color, s=200, alpha=0.8, edgecolors='black')
            ax.text(pos[node][0], pos[node][1], pos[node][2], node, fontsize=8)
        ax.set_title('3D Similarity Network')
        ax.set_axis_off()

# ============================================================================
# Interfacial Jamming Analyzer (with full PCA data export)
# ============================================================================
class InterfacialJammingAnalyzer
    
    Analyze signatures of interfacial jamming based on colloidal physics.

    FORMULAS & KEY REFERENCES
    ------------------------------------------------------------
    1. Jamming Index (J) J = G'  G'_max
       Trappe et al., Nature 411, 772-775 (2001); Nature Materials 7, 780-784 (2008).

    2. Fragility (F) = Payne Effect
       Payne, J. Appl. Polym. Sci. 6, 57-63 (1962).
       Hyun et al., Progress in Polymer Science 36, 1697-1753 (2011).

    3. Network Homogeneity (H) H = 1  (tanδ + 0.01)
       Larson, The Structure and Rheology of Complex Fluids, Oxford (1998).

    4. Critical Jamming Concentration (c)
       G' ~ A·(c - c)ᵝ
       Trappe et al., Nature Materials 7, 780-784 (2008).
       Zhang et al., PNAS 118, e2026481118 (2021).
    
    def __init__(self, parameters_df)
        self.df = parameters_df.copy()
        self.results = None
        self.critical_concentration = {}
        self.fit_params = {}

    def compute_jamming_indicators(self)
        print(n5. INTERFACIAL JAMMING ANALYSIS)
        print(-  40)

        # 1. Jamming index (normalized G')
        if 'G_prime_linear_Pa' in self.df.columns
            g_min = self.df['G_prime_linear_Pa'].min()
            g_max = self.df['G_prime_linear_Pa'].max()
            if g_max  g_min
                self.df['Jamming_index'] = (self.df['G_prime_linear_Pa'] - g_min)  (g_max - g_min)
            else
                self.df['Jamming_index'] = 0.5
        else
            self.df['Jamming_index'] = np.nan

        # 2. Fragility = Payne effect
        if 'Payne_effect' in self.df.columns
            self.df['Fragility'] = self.df['Payne_effect']
        else
            self.df['Fragility'] = np.nan

        # 3. Homogeneity = 1(tanδ + 0.01)
        if 'tan_delta_avg' in self.df.columns
            self.df['Homogeneity'] = 1  (self.df['tan_delta_avg'] + 0.01)
        else
            self.df['Homogeneity'] = np.nan

        # 4. Critical concentration estimation (power-law fit)
        self._estimate_jamming_concentrations()
        print(f  Estimated jamming concentrations 7S ~ {self.critical_concentration.get('7S', np.nan).2f}%, 11S ~ {self.critical_concentration.get('11S', np.nan).2f}%)

        # 5. Export initial table without PCA (will be updated after PCA)
        base_cols = ['Sample_ID', 'Protein', 'Concentration_%', 'G_prime_linear_Pa']
        optional_cols = ['Yield_stress_Pa', 'Viscosity_10s-1_Pas', 'Payne_effect',
                         'Jamming_index', 'Fragility', 'Homogeneity', 'tan_delta_avg']
        export_cols = [col for col in base_cols + optional_cols if col in self.df.columns]
        self.results = self.df[export_cols].copy()
        self.results.to_csv('Jamming_Indicators_temp.csv', index=False, encoding='utf-8-sig')
        print(  ✓ Saved Jamming_Indicators_temp.csv (will be overwritten with PCA scores later))

    def _estimate_jamming_concentrations(self)
        Fit G' = A(c - xc)^b to estimate critical concentration xc.
        for protein in ['7S', '11S']
            sub = self.df[self.df['Protein'] == protein].dropna(subset=['Concentration_%', 'G_prime_linear_Pa'])
            if len(sub) = 4
                x = sub['Concentration_%'].values
                y = sub['G_prime_linear_Pa'].values
                try
                    def power_law(c, A, xc, b)
                        return A  (c - xc)b
                    x_min = np.min(x)
                    popt, _ = curve_fit(power_law, x, y, p0=[1, x_min - 0.5, 1], maxfev=5000)
                    self.critical_concentration[protein] = popt[1]
                    self.fit_params[protein] = {'A' popt[0], 'xc' popt[1], 'beta' popt[2]}
                    # Add fitted values to dataframe
                    self.df.loc[sub.index, f'{protein}_Gprime_fit'] = power_law(x, popt)
                except Exception as e
                    print(f    Fit failed for {protein} {e})
                    self.critical_concentration[protein] = np.nan
                    self.fit_params[protein] = {}
            else
                self.critical_concentration[protein] = np.nan
                self.fit_params[protein] = {}

    def perform_jamming_pca(self)
        Perform PCA on jamming-related parameters and re-export data with PCA coordinates.
        print(n   Performing PCA on jamming metrics...)
        jam_params = ['G_prime_linear_Pa', 'Yield_stress_Pa', 'Payne_effect', 'tan_delta_avg', 'Recovery_percentage']
        available = [p for p in jam_params if p in self.df.columns]
        if not available
            print(  Not enough parameters for jamming PCA.)
            return

        X = self.df[available].fillna(self.df[available].median()).values
        X_scaled = StandardScaler().fit_transform(X)
        pca = PCA(n_components=2)
        pca_coords = pca.fit_transform(X_scaled)
        self.df['Jam_PC1'] = pca_coords[, 0]
        self.df['Jam_PC2'] = pca_coords[, 1]
        print(f    PCA explained variance PC1={pca.explained_variance_ratio_[0]100.1f}%, PC2={pca.explained_variance_ratio_[1]100.1f}%)

        # 2D plot
        fig, ax = plt.subplots(figsize=(8,6))
        for protein in ['7S', '11S']
            sub = self.df[self.df['Protein'] == protein]
            ax.scatter(sub['Jam_PC1'], sub['Jam_PC2'], color=PROTEIN_COLORS[protein],
                      s=100, alpha=0.8, edgecolors='black', label=protein)
            for _, row in sub.iterrows()
                ax.annotate(f{row['Concentration_%']}%, (row['Jam_PC1'], row['Jam_PC2']),
                            xytext=(5,5), textcoords='offset points', fontsize=8)
        ax.set_xlabel(f'Jam-PC1 ({pca.explained_variance_ratio_[0]100.1f}%)')
        ax.set_ylabel(f'Jam-PC2 ({pca.explained_variance_ratio_[1]100.1f}%)')
        ax.set_title('PCA of Jamming Metrics')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('Jamming_PCA.png', dpi=300)
        print(  ✓ Saved Jamming_PCA.png)

        # Re-export full jamming data with PCA coordinates
        self._export_full_jamming_data()

    def _export_full_jamming_data(self)
        Export final CSV files containing all jamming indicators + PCA scores.
        # Columns for Jamming_Raw_Data.csv (detailed)
        detailed_cols = ['Sample_ID', 'Protein', 'Concentration_%', 'G_prime_linear_Pa',
                         'Jam_PC1', 'Jam_PC2', 'Jamming_index', 'Fragility', 'Homogeneity',
                         'Payne_effect', 'tan_delta_avg', 'Yield_stress_Pa', 'Viscosity_10s-1_Pas',
                         'Recovery_percentage', 'Linear_limit_strain_%', 'Strain_stiffening_ratio']
        # Add fitted curve columns if exist
        fit_cols = [c for c in self.df.columns if 'Gprime_fit' in c]
        # Add critical concentration and fit parameters
        extra_cols = ['Critical_conc_7S_%', 'Critical_conc_11S_%']
        for protein in ['7S', '11S']
            if protein in self.fit_params and self.fit_params[protein]
                extra_cols += [f'{protein}_Fit_A', f'{protein}_Fit_xc', f'{protein}_Fit_beta']
        # Also add these to the dataframe as constant columns
        self.df['Critical_conc_7S_%'] = self.critical_concentration.get('7S', np.nan)
        self.df['Critical_conc_11S_%'] = self.critical_concentration.get('11S', np.nan)
        for protein in ['7S', '11S']
            if protein in self.fit_params
                self.df[f'{protein}_Fit_A'] = self.fit_params[protein].get('A', np.nan)
                self.df[f'{protein}_Fit_xc'] = self.fit_params[protein].get('xc', np.nan)
                self.df[f'{protein}_Fit_beta'] = self.fit_params[protein].get('beta', np.nan)

        all_cols = detailed_cols + fit_cols + extra_cols
        all_cols = [c for c in all_cols if c in self.df.columns]
        export_df = self.df[all_cols].copy()
        export_df.to_csv('Jamming_Raw_Data.csv', index=False, encoding='utf-8-sig')
        print(  ✓ Saved Jamming_Raw_Data.csv (includes Jam_PC1, Jam_PC2))

        # Simple version (Jamming_Indicators.csv)
        simple_cols = ['Sample_ID', 'Protein', 'Concentration_%', 'G_prime_linear_Pa',
                       'Jam_PC1', 'Jam_PC2', 'Jamming_index', 'Fragility', 'Homogeneity',
                       'Payne_effect', 'tan_delta_avg', 'Yield_stress_Pa', 'Viscosity_10s-1_Pas',
                       'Recovery_percentage']
        simple_cols = [c for c in simple_cols if c in self.df.columns]
        simple_df = self.df[simple_cols].copy()
        simple_df.to_csv('Jamming_Indicators.csv', index=False, encoding='utf-8-sig')
        print(  ✓ Saved Jamming_Indicators.csv (now includes Jam_PC1, Jam_PC2))

    def plot_jamming_diagrams(self)
        Plot jamming transition diagram (G' vs concentration) with power-law fit.
        fig, ax = plt.subplots(figsize=(8,6))
        for protein in ['7S', '11S']
            sub = self.df[self.df['Protein'] == protein].dropna(subset=['Concentration_%', 'G_prime_linear_Pa'])
            if len(sub)  0
                ax.plot(sub['Concentration_%'], sub['G_prime_linear_Pa'], 'o-',
                        color=PROTEIN_COLORS[protein], label=protein, markersize=10, linewidth=2)
                # Plot fitted curve if available
                fit_col = f'{protein}_Gprime_fit'
                if fit_col in self.df.columns and not self.df[fit_col].isna().all()
                    # Sort by concentration for smooth line
                    sub_fit = sub.sort_values('Concentration_%')
                    ax.plot(sub_fit['Concentration_%'], sub_fit[fit_col], '--',
                            color=PROTEIN_COLORS[protein], alpha=0.7, linewidth=1.5)
                # Mark critical point
                xc = self.critical_concentration.get(protein, np.nan)
                if not np.isnan(xc) and xc  min(sub['Concentration_%']) and xc  max(sub['Concentration_%'])
                    from scipy.interpolate import interp1d
                    f = interp1d(sub['Concentration_%'], sub['G_prime_linear_Pa'], kind='linear')
                    yc = f(xc)
                    ax.plot(xc, yc, 's', color=PROTEIN_COLORS[protein], markersize=12,
                            markerfacecolor='white', markeredgewidth=2)
                    ax.annotate(f'c ≈ {xc.1f}%', (xc, yc), xytext=(xc+0.3, yc1.1),
                                arrowprops=dict(arrowstyle='-', color='gray'))
        ax.set_xlabel('Protein Concentration (%)')
        ax.set_ylabel('Storage Modulus G' (Pa)')
        ax.set_title('Jamming Transition Diagram with Power-Law Fit')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('Jamming_Diagram.png', dpi=300)
        print(  ✓ Saved Jamming_Diagram.png)

# ============================================================================
# Report generation
# ============================================================================
def generate_report(df)
    with open('Analysis_Report_Protein.txt', 'w', encoding='utf-8') as f
        f.write(=  80 + n)
        f.write(RHEOLOGICAL ANALYSIS REPORT – SOY PROTEIN (7S & 11S)n)
        f.write(=  80 + nn)
        f.write(fAnalysis date {pd.Timestamp.now().strftime('%Y-%m-%d %H%M%S')}n)
        f.write(fNumber of samples {len(df)}n)
        f.write(fProtein types {df['Protein'].unique().tolist()}n)
        f.write(fConcentrations {sorted(df['Concentration_%'].unique())} %nn)
        if 'RPI' in df.columns
            best = df.loc[df['RPI'].idxmax(), 'Sample_ID']
            worst = df.loc[df['RPI'].idxmin(), 'Sample_ID']
            f.write(fBest performing sample {best} (RPI = {df['RPI'].max().3f})n)
            f.write(fWorst performing sample {worst} (RPI = {df['RPI'].min().3f})nn)
            f.write(Average RPI by proteinn)
            for prot in df['Protein'].unique()
                avg_rpi = df[df['Protein'] == prot]['RPI'].mean()
                f.write(f  {prot} {avg_rpi.3f}n)
        f.write(nCLUSTER ASSIGNMENTSn)
        f.write(-  40 + n)
        if 'Cluster' in df.columns
            for cl in sorted(df['Cluster'].unique())
                members = df[df['Cluster'] == cl]['Sample_ID'].tolist()
                f.write(fCluster {cl} {', '.join(members)}n)
        f.write(nRHEOLOGICAL PARAMETER SUMMARY (mean ± std)n)
        f.write(-  40 + n)
        for prot in df['Protein'].unique()
            sub = df[df['Protein'] == prot]
            f.write(fn{prot}n)
            for col in ['G_prime_linear_Pa', 'Yield_stress_Pa', 'Viscosity_10s-1_Pas', 'Recovery_percentage']
                if col in sub.columns
                    mean_val = sub[col].mean()
                    std_val = sub[col].std()
                    f.write(f  {col} {mean_val.2f} ± {std_val.2f}n)
        f.write(n + =  80 + n)
    print(  ✓ Saved Analysis_Report_Protein.txt)

# ============================================================================
# Export raw data for 2D3D plotting (general network PCA)
# ============================================================================
def export_plotting_data(df)
    print(n6. EXPORTING RAW DATA FOR 2D AND 3D PLOTTING)
    print(-  40)
    output = df.copy()
    if 'PC1_3D' not in output.columns
        numeric_cols = output.select_dtypes(include=[np.number]).columns.tolist()
        exclude = ['Concentration_%', 'RPI', 'Cluster', 'PC1', 'PC2', 'tSNE1', 'tSNE2']
        numeric_cols = [c for c in numeric_cols if c not in exclude]
        if len(numeric_cols) = 3
            X = output[numeric_cols].fillna(output[numeric_cols].median()).values
            X_scaled = StandardScaler().fit_transform(X)
            pca3 = PCA(n_components=3)
            pca3_coords = pca3.fit_transform(X_scaled)
            output['PC1_3D'] = pca3_coords[, 0]
            output['PC2_3D'] = pca3_coords[, 1]
            output['PC3_3D'] = pca3_coords[, 2]
            print(f  ✓ Computed 3D PCA coordinates (explained variance {pca3.explained_variance_ratio_}))
        else
            print(  Warning Not enough numeric parameters to compute 3D PCA.)
    output.to_csv('Rheology_Data_for_2D_3D_plots.csv', index=False, encoding='utf-8-sig')
    print(  ✓ Exported Rheology_Data_for_2D_3D_plots.csv)

# ============================================================================
# MAIN PIPELINE
# ============================================================================
def main()
    # 1. Load data
    loader = ProteinRheologyLoader(data_path='FDatarheology_data')
    data_dict = loader.load_all_data()
    if not data_dict
        print(No data loaded. Exiting.)
        return

    # 2. Extract parameters
    extractor = RheologyParameterExtractor(data_dict, loader.sample_info)
    param_df = extractor.extract_all_parameters()
    param_df.to_csv('Protein_Rheology_Parameters.csv', index=False)
    print(n✓ Saved Protein_Rheology_Parameters.csv)

    # 3. 2D mapping (adds PC1, PC2, tSNE1, tSNE2, Cluster)
    mapper = Network2DMapper(param_df)
    param_df = mapper.create_2d_maps()

    # 4. 3D visualizations
    visualizer = Network3DVisualizer(param_df)
    visualizer.create_3d_visualizations()

    # 5. Interfacial Jamming Analysis (compute indicators, then PCA, then export full data)
    jammer = InterfacialJammingAnalyzer(param_df)
    jammer.compute_jamming_indicators()
    jammer.plot_jamming_diagrams()
    jammer.perform_jamming_pca()   # This will also re-export CSV with PCA coordinates

    # 6. Export raw data for custom 2D3D network plotting (general)
    export_plotting_data(param_df)

    # 7. Generate report
    generate_report(param_df)

    # Print top performers
    print(nTop 3 samples by RPI)
    if 'RPI' in param_df.columns
        top3 = param_df.sort_values('RPI', ascending=False).head(3)
        for i, (_, row) in enumerate(top3.iterrows(), 1)
            print(f  {i}. {row['Sample_ID']} (RPI = {row['RPI'].3f}))

    print(n + =  80)
    print(ANALYSIS COMPLETE)
    print(=  80)

if __name__ == __main__
    main()