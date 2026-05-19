#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flow Rheology Comprehensive Analysis for Nature Communications Publication

This script performs comprehensive rheological analysis of 12 different salts,
incorporating 2D/3D visualization and molecular property correlations.
Exports raw data for external plotting.

Author: [Your Name]
Institution: [Your Institution]
"""

# ============================================================================
# IMPORTS
# ============================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.optimize import curve_fit
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
import warnings
warnings.filterwarnings('ignore')
from mpl_toolkits.mplot3d import Axes3D
import os
import traceback

# ============================================================================
# GLOBAL SETTINGS (Nature Communications Style)
# ============================================================================
plt.style.use('default')
plt.rcParams.update({
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
    'figure.autolayout': True,
    'figure.constrained_layout.use': False,
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.5,
    'lines.markersize': 6,
    'grid.linewidth': 0.5,
    'grid.alpha': 0.3,
    'font.family': ['Arial', 'DejaVu Sans'],
    'mathtext.fontset': 'stix',
})

print("=" * 80)
print("COMPREHENSIVE RHEOLOGICAL ANALYSIS")
print("Nature Communications Submission Format")
print("Exports raw data for external plotting")
print("=" * 80)

# ============================================================================
# MAIN ANALYSIS CLASS
# ============================================================================
class NatureCommsRheologyAnalyzer:
    """多sheet数据流变学分析器"""
    
    def __init__(self, data_path='F:/Data/rheology_data'):
        self.data_path = data_path
        self.all_data = {}
        self.salt_names = []
        self.combined_data = {}
        self.parameters_df = pd.DataFrame()
        self.molecular_properties = {}
        
        # 12种常见盐类名称
        self.standard_salts = [
            'NaOH', 'NaCl', 'KCl', 'CaCl2', 'MgCl2', 'AlCl3',
            'Na2CO3', 'NaHCO3', 'Na2SO4', 'K2SO4', 'MgSO4', 'CaSO4'
        ]
    
    def load_all_sheets(self):
        """加载所有CSV文件的所有sheet"""
        print("\n[STEP 1] Loading rheological data...")
        
        # 文件列表 - 添加分子性质文件
        files = {
            'laos': 'laos_data.csv',
            'freq': 'frequency_sweep.csv',
            'creep': 'creep_recovery.csv',
            'thix': 'thixotropy_loop.csv',
            'shear': 'steady_shear.csv',
            'molecular': 'molecular_descriptors.csv'  # 分子性质文件
        }
        
        # 首先尝试读取分子性质文件
        self._load_molecular_descriptors()
        
        # 然后尝试读取流变数据文件
        try:
            test_file = f'{self.data_path}/{files["laos"]}'
            
            # 尝试读取Excel
            try:
                xls = pd.ExcelFile(test_file)
                sheet_names = xls.sheet_names
                print(f"   Found {len(sheet_names)} sheets: {sheet_names}")
                self.salt_names = sheet_names
            except:
                # 如果不是Excel，可能是多个CSV文件
                print("   Trying to process as multiple CSV files...")
                self._load_as_multiple_csvs()
                return True
                
        except Exception as e:
            print(f"   File read error: {e}")
            return False
        
        # 读取所有流变数据文件的每个sheet
        for data_type, filename in files.items():
            if data_type == 'molecular':  # 分子性质已单独处理
                continue
                
            file_path = f'{self.data_path}/{filename}'
            
            if not os.path.exists(file_path):
                print(f"   ⚠ File not found: {filename}")
                continue
                
            print(f"\n   Processing: {filename}")
            self.all_data[data_type] = {}
            
            try:
                # 读取Excel文件的所有sheet
                xls = pd.ExcelFile(file_path)
                
                for sheet_name in xls.sheet_names:
                    try:
                        df = pd.read_excel(xls, sheet_name=sheet_name)
                        self.all_data[data_type][sheet_name] = df
                        print(f"    ✓ {sheet_name}: {df.shape[0]} rows, {df.shape[1]} cols")
                    except Exception as e:
                        print(f"    ✗ {sheet_name}: Failed - {e}")
                        
            except Exception as e:
                print(f"    Read failed: {e}")
                # 尝试作为普通CSV读取
                try:
                    df = pd.read_csv(file_path)
                    self.all_data[data_type]['Sheet1'] = df
                    print(f"    ✓ Read as CSV: {df.shape[0]} rows, {df.shape[1]} cols")
                except:
                    print(f"    ✗ Cannot read file")
        
        print(f"\n[RESULT] Data loaded, {len(self.all_data)} file types")
        return True
    
    def _load_molecular_descriptors(self):
        """加载分子性质描述符文件"""
        print("\n[STEP 1a] Loading molecular descriptors...")
        
        molecular_file = f'{self.data_path}/molecular_descriptors.csv'
        
        if not os.path.exists(molecular_file):
            print(f"   ⚠ Molecular descriptors file not found: {molecular_file}")
            print("   Will use default molecular properties")
            return
        
        try:
            # 尝试读取CSV文件
            mol_df = pd.read_csv(molecular_file)
            print(f"   ✓ Molecular descriptors loaded: {mol_df.shape[0]} salts, {mol_df.shape[1]} descriptors")
            
            # 显示列名以便调试
            print(f"   Columns in molecular descriptors: {list(mol_df.columns)}")
            
            # 确定盐类名称列 - 尝试多种可能的列名
            salt_column = None
            for col in ['Salt', 'salt', '盐类', '盐名称', '名称', 'Name', 'name']:
                if col in mol_df.columns:
                    salt_column = col
                    break
            
            if salt_column is None:
                print("   ⚠ Could not find salt name column in molecular descriptors")
                # 假设第一列是盐类名称
                salt_column = mol_df.columns[0]
                print(f"   Using first column as salt names: {salt_column}")
            
            # 处理每一行的分子性质
            for idx, row in mol_df.iterrows():
                salt_name = str(row[salt_column]).strip()
                
                # 跳过空值
                if pd.isna(salt_name) or salt_name == '':
                    continue
                
                # 创建分子性质字典
                mol_props = {}
                
                # 尝试提取常见的分子性质
                properties_to_extract = {
                    'LogP': ['LogP', 'logP', 'logp', 'Log P', 'log P', '辛醇水分配系数'],
                    'TPSA': ['TPSA', 'tpsa', 'TPSA_1', 'TPSA_A', '拓扑极性表面积', '极性表面积'],
                    'MW': ['MW', 'mw', 'MolecularWeight', 'molecular_weight', '分子量', '分子重量'],
                    'HBD': ['HBD', 'hbd', 'HydrogenBondDonor', '氢键供体', '供氢原子'],
                    'HBA': ['HBA', 'hba', 'HydrogenBondAcceptor', '氢键受体', '受氢原子'],
                    'RotBonds': ['RotBonds', 'rotatable_bonds', 'RotatableBonds', '可旋转键', '旋转键'],
                    'Charge': ['Charge', 'charge', '电荷', '离子电荷'],
                    'IonicRadius': ['IonicRadius', 'ionic_radius', '离子半径', '半径'],
                    'HydrationEnergy': ['HydrationEnergy', 'hydration_energy', '水合能', '水合能量']
                }
                
                for prop_name, possible_names in properties_to_extract.items():
                    for name in possible_names:
                        if name in mol_df.columns:
                            value = row[name]
                            if not pd.isna(value):
                                try:
                                    mol_props[prop_name] = float(value)
                                except:
                                    mol_props[prop_name] = value
                            break
                
                # 如果有提取到分子性质，则保存
                if mol_props:
                    self.molecular_properties[salt_name] = mol_props
                    print(f"    ✓ {salt_name}: {list(mol_props.keys())}")
                else:
                    print(f"    ⚠ {salt_name}: No molecular properties extracted")
            
            print(f"   Total molecular properties loaded for {len(self.molecular_properties)} salts")
            
        except Exception as e:
            print(f"   ✗ Error loading molecular descriptors: {e}")
            print(traceback.format_exc())
    
    def _load_as_multiple_csvs(self):
        """如果文件是多个CSV的处理方法"""
        print("   Trying to find individual CSV files for 12 salts...")
        
        for salt in self.standard_salts:
            for data_type in ['laos', 'freq', 'shear']:
                filename = f'{salt}_{data_type}.csv'
                file_path = f'{self.data_path}/{filename}'
                
                if os.path.exists(file_path):
                    try:
                        df = pd.read_csv(file_path)
                        if data_type not in self.all_data:
                            self.all_data[data_type] = {}
                        self.all_data[data_type][salt] = df
                        print(f"    ✓ {salt}_{data_type}: {df.shape}")
                    except Exception as e:
                        print(f"    ✗ {salt}_{data_type}: {e}")
        
        self.salt_names = list(set([salt for data_type in self.all_data.values() 
                                   for salt in data_type.keys()]))
        print(f"\n   Found {len(self.salt_names)} salts: {self.salt_names}")
        
        # 加载分子性质（如果还未加载）
        if not self.molecular_properties:
            self._load_molecular_descriptors()
    
    def combine_data_by_salt(self):
        """按盐类合并数据"""
        print("\n[STEP 2] Combining data by salt...")
        
        if not self.all_data:
            print("   No data available")
            return
        
        self.combined_data = {}
        
        # 首先从流变数据中收集所有盐类名称
        all_salts_from_rheology = set()
        for data_type, salt_dict in self.all_data.items():
            for salt_name in salt_dict.keys():
                all_salts_from_rheology.add(salt_name)
        
        # 添加分子性质中的盐类
        all_salts_from_rheology.update(self.molecular_properties.keys())
        
        self.salt_names = list(all_salts_from_rheology)
        print(f"   Total salts identified: {len(self.salt_names)}")
        
        # 初始化数据结构
        for salt in self.salt_names:
            self.combined_data[salt] = {}
        
        # 合并流变数据
        for data_type, salt_dict in self.all_data.items():
            for salt_name, df in salt_dict.items():
                if salt_name in self.combined_data:
                    self.combined_data[salt_name][data_type] = df
                    print(f"  {salt_name}: Added {data_type} data")
        
        # 添加分子性质数据
        for salt_name, mol_props in self.molecular_properties.items():
            if salt_name in self.combined_data:
                self.combined_data[salt_name]['molecular'] = pd.DataFrame([mol_props])
                print(f"  {salt_name}: Added molecular properties")
        
        print(f"\n[RESULT] Combined data for {len(self.combined_data)} salts")
    
    def extract_rheological_parameters(self):
        """为每种盐类提取流变学参数"""
        print("\n[STEP 3] Extracting rheological parameters...")
        
        salt_parameters = []
        
        for salt_name, data_dict in self.combined_data.items():
            params = {'盐类': salt_name}
            
            # 1. LAOS参数
            if 'laos' in data_dict:
                laos_data = data_dict['laos']
                
                if 'storage_modulus(Pa)' in laos_data.columns:
                    G_prime_initial = laos_data['storage_modulus(Pa)'].iloc[0]
                    params['G_prime_linear'] = G_prime_initial
                    
                    if len(laos_data) > 1:
                        G_prime_min = laos_data['storage_modulus(Pa)'].min()
                        if G_prime_initial > 0:
                            params['Payne_effect'] = (G_prime_initial - G_prime_min) / G_prime_initial
                    
                    G_prime_max = laos_data['storage_modulus(Pa)'].max()
                    if G_prime_initial > 0:
                        params['Strain_stiffening'] = G_prime_max / G_prime_initial
            
            # 2. 频率扫描参数
            if 'freq' in data_dict:
                freq_data = data_dict['freq']
                
                if 'storage_modulus(Pa)' in freq_data.columns and 'loss_modulus(Pa)' in freq_data.columns:
                    params['G_prime_avg'] = freq_data['storage_modulus(Pa)'].mean()
                    params['G_double_prime_avg'] = freq_data['loss_modulus(Pa)'].mean()
                    
                    if params['G_prime_avg'] > 0:
                        params['tan_delta_avg'] = params['G_double_prime_avg'] / params['G_prime_avg']
            
            # 3. 稳态剪切参数
            if 'shear' in data_dict:
                shear_data = data_dict['shear']
                
                if 'viscosity(Pa.s)' in shear_data.columns:
                    params['avg_viscosity'] = shear_data['viscosity(Pa.s)'].mean()
                
                if 'stress(Pa)' in shear_data.columns:
                    params['yield_stress'] = shear_data['stress(Pa)'].iloc[0]
                    
                    if 'shear_rate(1/s)' in shear_data.columns:
                        try:
                            rates = shear_data['shear_rate(1/s)'].values
                            stresses = shear_data['stress(Pa)'].values
                            
                            valid = (rates > 0) & (stresses > 0)
                            if np.sum(valid) > 2:
                                log_rates = np.log10(rates[valid])
                                log_stresses = np.log10(stresses[valid])
                                
                                coeffs = np.polyfit(log_rates, log_stresses, 1)
                                params['flow_index'] = coeffs[0]
                                params['consistency_index'] = 10**coeffs[1]
                        except:
                            pass
            
            # 4. 蠕变恢复参数
            if 'creep' in data_dict:
                creep_data = data_dict['creep']
                if 'strain(%)' in creep_data.columns:
                    max_strain = creep_data['strain(%)'].max()
                    final_strain = creep_data['strain(%)'].iloc[-1]
                    
                    if max_strain > 0:
                        params['recovery_percentage'] = (1 - final_strain / max_strain) * 100
            
            # 5. 添加分子性质
            if salt_name in self.molecular_properties:
                params.update(self.molecular_properties[salt_name])
            
            salt_parameters.append(params)
        
        # 创建参数DataFrame
        self.parameters_df = pd.DataFrame(salt_parameters)
        
        # 计算综合性能指数
        self._calculate_performance_index()
        
        # 计算复合模量
        if 'G_prime_avg' in self.parameters_df.columns and 'G_double_prime_avg' in self.parameters_df.columns:
            self.parameters_df['G_star'] = np.sqrt(
                self.parameters_df['G_prime_avg']**2 + 
                self.parameters_df['G_double_prime_avg']**2
            )
        
        print(f"\n[RESULT] Extracted {len(self.parameters_df.columns)} parameters for {len(self.parameters_df)} salts")
        print(f"   Parameters: {list(self.parameters_df.columns)}")
        
        return self.parameters_df
    
    def _calculate_performance_index(self):
        """计算综合性能指数"""
        key_params = []
        for param in ['G_prime_linear', 'avg_viscosity', 'recovery_percentage', 'Strain_stiffening']:
            if param in self.parameters_df.columns:
                key_params.append(param)
        
        if not key_params:
            return
        
        def normalize(series):
            min_val = series.min()
            max_val = series.max()
            if max_val > min_val:
                return (series - min_val) / (max_val - min_val)
            else:
                return pd.Series([0.5] * len(series))
        
        for param in key_params:
            self.parameters_df[f'{param}_norm'] = normalize(self.parameters_df[param].fillna(
                self.parameters_df[param].median()))
        
        rpi_values = []
        for idx, row in self.parameters_df.iterrows():
            rpi = 0
            count = 0
            for param in key_params:
                norm_col = f'{param}_norm'
                if norm_col in row.index and not pd.isna(row[norm_col]):
                    rpi += row[norm_col]
                    count += 1
            if count > 0:
                rpi_values.append(rpi / count)
            else:
                rpi_values.append(0.5)
        
        self.parameters_df['RPI'] = rpi_values
    
    def analyze_salt_effects(self):
        """分析盐类效应"""
        print("\n[STEP 4] Analyzing salt effects...")
        
        if self.parameters_df.empty:
            print("   No parameter data available")
            return
        
        # 创建分析图
        fig = plt.figure(figsize=(15, 10))
        
        # 1. RPI排名
        ax1 = plt.subplot(2, 3, 1)
        sorted_df = self.parameters_df.sort_values('RPI', ascending=True)
        bars = ax1.barh(range(len(sorted_df)), sorted_df['RPI'], 
                       color=plt.cm.viridis(np.linspace(0, 1, len(sorted_df))))
        ax1.set_yticks(range(len(sorted_df)))
        ax1.set_yticklabels(sorted_df['盐类'])
        ax1.set_xlabel('Rheological Performance Index (RPI)')
        ax1.set_title('Salt Performance Ranking')
        ax1.grid(True, alpha=0.3, axis='x')
        
        # 2. 储能模量对比
        ax2 = plt.subplot(2, 3, 2)
        if 'G_prime_linear' in self.parameters_df.columns:
            sorted_gprime = self.parameters_df.sort_values('G_prime_linear', ascending=True)
            bars = ax2.barh(range(len(sorted_gprime)), sorted_gprime['G_prime_linear'],
                           color=plt.cm.plasma(np.linspace(0, 1, len(sorted_gprime))))
            ax2.set_yticks(range(len(sorted_gprime)))
            ax2.set_yticklabels(sorted_gprime['盐类'])
            ax2.set_xlabel('Storage Modulus G\' (Pa)')
            ax2.set_title('Elasticity Comparison')
            ax2.grid(True, alpha=0.3, axis='x')
        
        # 3. 粘度对比
        ax3 = plt.subplot(2, 3, 3)
        if 'avg_viscosity' in self.parameters_df.columns:
            sorted_visc = self.parameters_df.sort_values('avg_viscosity', ascending=True)
            bars = ax3.barh(range(len(sorted_visc)), sorted_visc['avg_viscosity'],
                           color=plt.cm.cool(np.linspace(0, 1, len(sorted_visc))))
            ax3.set_yticks(range(len(sorted_visc)))
            ax3.set_yticklabels(sorted_visc['盐类'])
            ax3.set_xlabel('Average Viscosity (Pa·s)')
            ax3.set_title('Viscosity Comparison')
            ax3.grid(True, alpha=0.3, axis='x')
        
        # 4. 二维性能矩阵
        ax4 = plt.subplot(2, 3, 4)
        if 'G_prime_linear' in self.parameters_df.columns and 'avg_viscosity' in self.parameters_df.columns:
            scatter = ax4.scatter(self.parameters_df['G_prime_linear'], 
                                 self.parameters_df['avg_viscosity'],
                                 c=self.parameters_df['RPI'], s=150, 
                                 cmap='viridis', alpha=0.8, edgecolors='black')
            
            for idx, row in self.parameters_df.iterrows():
                ax4.annotate(row['盐类'], 
                            xy=(row['G_prime_linear'], row['avg_viscosity']),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=9, fontweight='bold')
            
            ax4.set_xlabel('Storage Modulus G\' (Pa)')
            ax4.set_ylabel('Average Viscosity (Pa·s)')
            ax4.set_title('Performance Matrix (color: RPI)')
            ax4.grid(True, alpha=0.3)
            
            # 添加colorbar
            cbar_ax = fig.add_axes([0.15, 0.05, 0.2, 0.02])
            plt.colorbar(scatter, cax=cbar_ax, orientation='horizontal')
            cbar_ax.set_xlabel('RPI')
        
        # 5. LogP与RPI关系
        ax5 = plt.subplot(2, 3, 5)
        if 'LogP' in self.parameters_df.columns and 'RPI' in self.parameters_df.columns:
            ax5.scatter(self.parameters_df['LogP'], self.parameters_df['RPI'],
                       s=100, alpha=0.7, edgecolors='black')
            
            for idx, row in self.parameters_df.iterrows():
                ax5.annotate(row['盐类'], 
                           xy=(row['LogP'], row['RPI']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=9, fontweight='bold')
            
            ax5.set_xlabel('LogP')
            ax5.set_ylabel('RPI')
            ax5.set_title('LogP vs Rheological Performance')
            ax5.grid(True, alpha=0.3)
        
        # 6. TPSA与RPI关系
        ax6 = plt.subplot(2, 3, 6)
        if 'TPSA' in self.parameters_df.columns and 'RPI' in self.parameters_df.columns:
            ax6.scatter(self.parameters_df['TPSA'], self.parameters_df['RPI'],
                       s=100, alpha=0.7, edgecolors='black')
            
            for idx, row in self.parameters_df.iterrows():
                ax6.annotate(row['盐类'], 
                           xy=(row['TPSA'], row['RPI']),
                           xytext=(5, 5), textcoords='offset points',
                           fontsize=9, fontweight='bold')
            
            ax6.set_xlabel('TPSA')
            ax6.set_ylabel('RPI')
            ax6.set_title('TPSA vs Rheological Performance')
            ax6.grid(True, alpha=0.3)
        
        # 手动调整布局
        plt.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.12, 
                          wspace=0.3, hspace=0.4)
        
        plt.savefig('Figure_Salt_Comparison.png', dpi=300, bbox_inches='tight')
        print("   ✓ Figure saved: Figure_Salt_Comparison.png")
        
        return fig
    
    def create_2d_network_analysis(self):
        """创建二维网络分析 - 包含分子性质"""
        print("\n[STEP 5] Creating 2D network analysis...")
        
        if self.parameters_df.empty:
            print("   No parameter data available")
            return
        
        # 选择数值型特征 - 包括分子性质
        numeric_cols = self.parameters_df.select_dtypes(include=[np.number]).columns.tolist()
        exclude_cols = ['RPI']
        feature_cols = [c for c in numeric_cols if c not in exclude_cols]
        
        if len(feature_cols) < 2:
            print("   Insufficient features")
            return
        
        print(f"   Features for PCA: {len(feature_cols)} numerical columns")
        print(f"   Including molecular properties: {[c for c in feature_cols if c in ['LogP', 'TPSA', 'MW', 'HBD', 'HBA', 'RotBonds', 'Charge', 'IonicRadius', 'HydrationEnergy']]}")
        
        # 填充缺失值
        data_for_pca = self.parameters_df[feature_cols].fillna(
            self.parameters_df[feature_cols].median())
        
        # PCA降维
        pca = PCA(n_components=2)
        scaled_data = StandardScaler().fit_transform(data_for_pca)
        pca_result = pca.fit_transform(scaled_data)
        
        self.parameters_df['PC1'] = pca_result[:, 0]
        self.parameters_df['PC2'] = pca_result[:, 1]
        
        # 聚类分析
        kmeans = KMeans(n_clusters=min(4, len(self.parameters_df)), random_state=42, n_init=10)
        self.parameters_df['Cluster'] = kmeans.fit_predict(scaled_data)
        
        # 创建单独的图形，避免tight_layout问题
        # 图形1: PCA散点图（按RPI着色）
        fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        fig1.suptitle('PCA Analysis of Rheological and Molecular Parameters', fontsize=14, fontweight='bold')
        
        # 1. 按RPI着色
        scatter1 = ax1.scatter(self.parameters_df['PC1'], self.parameters_df['PC2'],
                              c=self.parameters_df['RPI'], s=150, cmap='viridis',
                              alpha=0.8, edgecolors='black')
        
        for idx, row in self.parameters_df.iterrows():
            ax1.annotate(row['盐类'], xy=(row['PC1'], row['PC2']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=9, fontweight='bold')
        
        ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        ax1.set_title('Colored by RPI')
        ax1.grid(True, alpha=0.3)
        
        # 创建独立的colorbar轴
        cbar_ax1 = fig1.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar1 = plt.colorbar(scatter1, cax=cbar_ax1)
        cbar1.set_label('RPI', rotation=270, labelpad=15)
        
        # 2. 按聚类着色
        scatter2 = ax2.scatter(self.parameters_df['PC1'], self.parameters_df['PC2'],
                              c=self.parameters_df['Cluster'], s=150, cmap='tab20',
                              alpha=0.8, edgecolors='black')
        
        for idx, row in self.parameters_df.iterrows():
            ax2.annotate(row['盐类'], xy=(row['PC1'], row['PC2']),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=9, fontweight='bold')
        
        ax2.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
        ax2.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
        ax2.set_title(f'Colored by Cluster (k={min(4, len(self.parameters_df))})')
        ax2.grid(True, alpha=0.3)
        
        # 为聚类图添加colorbar
        cbar_ax2 = fig1.add_axes([0.95, 0.15, 0.02, 0.7])
        cbar2 = plt.colorbar(scatter2, cax=cbar_ax2)
        cbar2.set_label('Cluster', rotation=270, labelpad=15)
        
        # 调整布局
        plt.subplots_adjust(left=0.08, right=0.9, top=0.9, bottom=0.1, wspace=0.25)
        
        plt.savefig('Figure_2D_PCA_Analysis.png', dpi=300, bbox_inches='tight')
        print("   ✓ Figure saved: Figure_2D_PCA_Analysis.png")
        
        # 图形2: 特征贡献度和聚类分析
        fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(14, 6))
        fig2.suptitle('Feature Importance and Cluster Analysis', fontsize=14, fontweight='bold')
        
        # 3. 特征贡献度
        feature_importance = pd.DataFrame({
            'feature': feature_cols,
            'importance': np.abs(pca.components_[0])
        }).sort_values('importance', ascending=False).head(10)
        
        bars = ax3.barh(range(len(feature_importance)), feature_importance['importance'],
                       color=plt.cm.coolwarm(np.linspace(0, 1, len(feature_importance))))
        ax3.set_yticks(range(len(feature_importance)))
        ax3.set_yticklabels(feature_importance['feature'])
        ax3.set_xlabel('Contribution to PC1')
        ax3.set_title('Top 10 Feature Contributions')
        ax3.grid(True, alpha=0.3, axis='x')
        
        # 4. 聚类性能分析
        if 'Cluster' in self.parameters_df.columns:
            cluster_stats = self.parameters_df.groupby('Cluster').agg({
                'RPI': 'mean',
                'G_prime_linear': 'mean' if 'G_prime_linear' in self.parameters_df.columns else None,
                'avg_viscosity': 'mean' if 'avg_viscosity' in self.parameters_df.columns else None
            }).dropna(axis=1, how='all')
            
            if not cluster_stats.empty:
                cluster_stats.plot(kind='bar', ax=ax4, alpha=0.7)
                ax4.set_xlabel('Cluster')
                ax4.set_ylabel('Mean Value')
                ax4.set_title('Cluster Performance Comparison')
                ax4.legend()
                ax4.grid(True, alpha=0.3)
        
        # 调整布局
        plt.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.1, wspace=0.3)
        
        plt.savefig('Figure_Feature_Cluster_Analysis.png', dpi=300, bbox_inches='tight')
        print("   ✓ Figure saved: Figure_Feature_Cluster_Analysis.png")
        
        # 显示聚类结果
        print("\n   Clustering Results:")
        for cluster in sorted(self.parameters_df['Cluster'].unique()):
            cluster_salts = self.parameters_df[self.parameters_df['Cluster'] == cluster]['盐类'].tolist()
            cluster_rpi = self.parameters_df[self.parameters_df['Cluster'] == cluster]['RPI'].mean()
            print(f"     Cluster {cluster}: {cluster_salts} (Avg RPI: {cluster_rpi:.3f})")
        
        return fig1, fig2
    
    def create_3d_visualization(self):
        """创建三维可视化 - 包含分子性质"""
        print("\n[STEP 6] Creating 3D visualization...")
        
        if self.parameters_df.empty:
            print("   No parameter data available")
            return
        
        # 选择三个关键参数 - 尝试包括分子性质
        key_params = []
        
        # 优先选择分子性质
        molecular_params = ['LogP', 'TPSA', 'MW', 'Charge', 'IonicRadius']
        rheology_params = ['G_prime_linear', 'avg_viscosity', 'tan_delta_avg', 'Payne_effect', 'flow_index']
        
        # 先添加分子性质
        for param in molecular_params:
            if param in self.parameters_df.columns and len(key_params) < 3:
                key_params.append(param)
        
        # 如果分子性质不足，添加流变参数
        for param in rheology_params:
            if param in self.parameters_df.columns and len(key_params) < 3:
                key_params.append(param)
        
        if len(key_params) < 3:
            print(f"   Insufficient key parameters (have {len(key_params)}, need 3)")
            # 使用前三个可用的数值列
            numeric_cols = self.parameters_df.select_dtypes(include=[np.number]).columns.tolist()
            if len(numeric_cols) >= 3:
                key_params = numeric_cols[:3]
            else:
                return
        
        print(f"   3D visualization parameters: {key_params}")
        
        # 创建单独的3D图形
        fig1 = plt.figure(figsize=(10, 8))
        ax1 = fig1.add_subplot(111, projection='3d')
        
        x = self.parameters_df[key_params[0]].fillna(0).values
        y = self.parameters_df[key_params[1]].fillna(0).values
        z = self.parameters_df[key_params[2]].fillna(0).values
        
        # 根据RPI着色
        scatter = ax1.scatter(x, y, z, c=self.parameters_df['RPI'], 
                             cmap='viridis', s=150, alpha=0.8, edgecolors='black')
        
        # 添加标签
        for i, row in self.parameters_df.iterrows():
            ax1.text(x[i], y[i], z[i], row['盐类'], fontsize=9, fontweight='bold')
        
        ax1.set_xlabel(key_params[0])
        ax1.set_ylabel(key_params[1])
        ax1.set_zlabel(key_params[2])
        ax1.set_title('3D Parameter Space (Colored by RPI)')
        
        # 添加colorbar到右侧
        cbar_ax = fig1.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = plt.colorbar(scatter, cax=cbar_ax)
        cbar.set_label('RPI', rotation=270, labelpad=15)
        
        plt.savefig('Figure_3D_Parameter_Space.png', dpi=300, bbox_inches='tight')
        print("   ✓ Figure saved: Figure_3D_Parameter_Space.png")
        
        # 创建平行坐标图
        fig2, ax2 = plt.subplots(figsize=(12, 8))
        
        parallel_params = ['RPI']
        
        # 添加流变参数
        for param in ['G_prime_linear', 'avg_viscosity', 'tan_delta_avg', 'Payne_effect']:
            if param in self.parameters_df.columns and len(parallel_params) < 5:
                parallel_params.append(param)
        
        # 添加分子性质
        for param in ['LogP', 'TPSA', 'MW', 'Charge']:
            if param in self.parameters_df.columns and len(parallel_params) < 8:
                parallel_params.append(param)
        
        if len(parallel_params) > 1:
            data_to_plot = self.parameters_df[parallel_params].copy()
            for col in data_to_plot.columns:
                data_to_plot[col] = data_to_plot[col].fillna(data_to_plot[col].median())
            
            scaler = MinMaxScaler()
            data_normalized = scaler.fit_transform(data_to_plot)
            data_normalized_df = pd.DataFrame(data_normalized, columns=parallel_params)
            
            colors = plt.cm.viridis(self.parameters_df['RPI'].values / self.parameters_df['RPI'].max())
            
            for i in range(len(self.parameters_df)):
                ax2.plot(range(len(parallel_params)), data_normalized_df.iloc[i], 
                        color=colors[i], alpha=0.7, linewidth=2)
            
            ax2.set_xticks(range(len(parallel_params)))
            ax2.set_xticklabels(parallel_params, rotation=45, ha='right')
            ax2.set_ylabel('Normalized Value')
            ax2.set_title('Parallel Coordinates: Rheological and Molecular Properties')
            ax2.grid(True, alpha=0.3)
            
            # 添加图例
            legend_elements = []
            top_salts = self.parameters_df.sort_values('RPI', ascending=False).head(5)
            for i, (_, row) in enumerate(top_salts.iterrows()):
                color_idx = list(self.parameters_df['盐类']).index(row['盐类'])
                legend_elements.append(plt.Line2D([0], [0], color=colors[color_idx], 
                                                 lw=2, label=f"{row['盐类']} (RPI: {row['RPI']:.2f})"))
            
            ax2.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1))
            
            # 调整布局给图例留空间
            plt.subplots_adjust(right=0.75)
        
        plt.savefig('Figure_Parallel_Coordinates.png', dpi=300, bbox_inches='tight')
        print("   ✓ Figure saved: Figure_Parallel_Coordinates.png")
        
        return fig1, fig2
    
    def create_molecular_correlation_analysis(self):
        """创建分子性质与流变性能的相关性分析"""
        print("\n[STEP 7] Creating molecular property - rheology correlation analysis...")
        
        if self.parameters_df.empty:
            print("   No parameter data available")
            return
        
        # 检查是否有分子性质数据
        molecular_props = ['LogP', 'TPSA', 'MW', 'HBD', 'HBA', 'RotBonds', 'Charge', 'IonicRadius', 'HydrationEnergy']
        available_molecular = [p for p in molecular_props if p in self.parameters_df.columns]
        
        if len(available_molecular) == 0:
            print("   No molecular property data available")
            return
        
        # 检查是否有流变参数
        rheology_params = ['RPI', 'G_prime_linear', 'avg_viscosity', 'tan_delta_avg', 'Payne_effect', 'flow_index', 'consistency_index']
        available_rheology = [p for p in rheology_params if p in self.parameters_df.columns]
        
        if len(available_rheology) == 0:
            print("   No rheological parameter data available")
            return
        
        print(f"   Available molecular properties: {available_molecular}")
        print(f"   Available rheological parameters: {available_rheology}")
        
        # 创建相关性图形
        fig = plt.figure(figsize=(14, 10))
        fig.suptitle('Molecular Property - Rheology Correlations', fontsize=14, fontweight='bold')
        
        # 1. 分子性质与RPI的相关性
        ax1 = plt.subplot(2, 2, 1)
        if 'RPI' in self.parameters_df.columns:
            corr_values = []
            p_values = []
            
            for mol_prop in available_molecular:
                clean_data = self.parameters_df[[mol_prop, 'RPI']].dropna()
                if len(clean_data) >= 3:
                    corr, p_val = stats.pearsonr(clean_data[mol_prop], clean_data['RPI'])
                    corr_values.append(corr)
                    p_values.append(p_val)
                else:
                    corr_values.append(0)
                    p_values.append(1)
            
            # 根据p值确定颜色
            colors = []
            for p_val in p_values:
                if p_val < 0.001:
                    colors.append('#2ecc71')  # 绿色：极显著
                elif p_val < 0.01:
                    colors.append('#3498db')  # 蓝色：非常显著
                elif p_val < 0.05:
                    colors.append('#f39c12')  # 橙色：显著
                else:
                    colors.append('#e74c3c')  # 红色：不显著
            
            bars = ax1.bar(range(len(available_molecular)), corr_values, color=colors)
            
            ax1.set_xticks(range(len(available_molecular)))
            ax1.set_xticklabels(available_molecular, rotation=45)
            ax1.set_ylabel('Correlation Coefficient')
            ax1.set_title('Molecular Properties vs RPI')
            ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
            ax1.grid(True, alpha=0.3, axis='y')
            
            # 添加数值和显著性标记
            for i, (corr, p_val) in enumerate(zip(corr_values, p_values)):
                significance = ''
                if p_val < 0.001:
                    significance = '***'
                elif p_val < 0.01:
                    significance = '**'
                elif p_val < 0.05:
                    significance = '*'
                
                ax1.text(i, corr + (0.02 if corr >= 0 else -0.05), 
                        f'{corr:.3f}{significance}', 
                        ha='center', va='bottom' if corr >= 0 else 'top', 
                        fontsize=8, fontweight='bold')
        
        # 2. 相关性矩阵热图
        ax2 = plt.subplot(2, 2, 2)
        
        # 组合用于相关性分析的列
        corr_columns = available_molecular + available_rheology
        corr_data = self.parameters_df[corr_columns].fillna(0)
        corr_matrix = corr_data.corr()
        
        # 提取分子性质和流变参数之间的相关性
        molecular_rheo_corr = corr_matrix.loc[available_molecular, available_rheology]
        
        im = ax2.imshow(molecular_rheo_corr.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        
        # 添加数值
        for i in range(len(molecular_rheo_corr)):
            for j in range(len(molecular_rheo_corr.columns)):
                value = molecular_rheo_corr.iloc[i, j]
                if not np.isnan(value):
                    color = 'white' if abs(value) > 0.5 else 'black'
                    ax2.text(j, i, f'{value:.2f}', ha='center', va='center', 
                            color=color, fontsize=8)
        
        ax2.set_xticks(range(len(molecular_rheo_corr.columns)))
        ax2.set_yticks(range(len(molecular_rheo_corr)))
        ax2.set_xticklabels(molecular_rheo_corr.columns, rotation=45, ha='right', fontsize=8)
        ax2.set_yticklabels(molecular_rheo_corr.index, fontsize=8)
        ax2.set_title('Molecular-Rheology Correlation Matrix')
        
        # 添加colorbar
        cbar_ax = fig.add_axes([0.92, 0.55, 0.02, 0.3])
        plt.colorbar(im, cax=cbar_ax)
        cbar_ax.set_ylabel('Correlation', rotation=270, labelpad=15)
        
        # 3. LogP与G'的关系
        ax3 = plt.subplot(2, 2, 3)
        if 'LogP' in self.parameters_df.columns and 'G_prime_linear' in self.parameters_df.columns:
            clean_data = self.parameters_df[['LogP', 'G_prime_linear']].dropna()
            if len(clean_data) >= 3:
                ax3.scatter(clean_data['LogP'], clean_data['G_prime_linear'], 
                           s=100, alpha=0.7, edgecolors='black')
                
                # 添加回归线
                try:
                    z = np.polyfit(clean_data['LogP'], clean_data['G_prime_linear'], 1)
                    p = np.poly1d(z)
                    ax3.plot(clean_data['LogP'], p(clean_data['LogP']), 
                            "r--", alpha=0.8, label=f'y={z[0]:.2f}x+{z[1]:.2f}')
                    
                    # 计算R²
                    residuals = clean_data['G_prime_linear'] - p(clean_data['LogP'])
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((clean_data['G_prime_linear'] - clean_data['G_prime_linear'].mean())**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    
                    ax3.legend(title=f'R² = {r_squared:.3f}')
                except:
                    pass
                
                # 添加标签
                for i, row in clean_data.iterrows():
                    salt_name = self.parameters_df.loc[i, '盐类']
                    ax3.text(row['LogP'], row['G_prime_linear'], 
                            salt_name, fontsize=8, alpha=0.7)
                
                ax3.set_xlabel('LogP')
                ax3.set_ylabel('G\' (Pa)')
                ax3.set_title('LogP vs Storage Modulus')
                ax3.grid(True, alpha=0.3)
        
        # 4. TPSA与粘度的关系
        ax4 = plt.subplot(2, 2, 4)
        if 'TPSA' in self.parameters_df.columns and 'avg_viscosity' in self.parameters_df.columns:
            clean_data = self.parameters_df[['TPSA', 'avg_viscosity']].dropna()
            if len(clean_data) >= 3:
                ax4.scatter(clean_data['TPSA'], clean_data['avg_viscosity'], 
                           s=100, alpha=0.7, edgecolors='black')
                
                # 添加回归线
                try:
                    z = np.polyfit(clean_data['TPSA'], clean_data['avg_viscosity'], 1)
                    p = np.poly1d(z)
                    ax4.plot(clean_data['TPSA'], p(clean_data['TPSA']), 
                            "r--", alpha=0.8, label=f'y={z[0]:.2f}x+{z[1]:.2f}')
                    
                    # 计算R²
                    residuals = clean_data['avg_viscosity'] - p(clean_data['TPSA'])
                    ss_res = np.sum(residuals**2)
                    ss_tot = np.sum((clean_data['avg_viscosity'] - clean_data['avg_viscosity'].mean())**2)
                    r_squared = 1 - (ss_res / ss_tot)
                    
                    ax4.legend(title=f'R² = {r_squared:.3f}')
                except:
                    pass
                
                # 添加标签
                for i, row in clean_data.iterrows():
                    salt_name = self.parameters_df.loc[i, '盐类']
                    ax4.text(row['TPSA'], row['avg_viscosity'], 
                            salt_name, fontsize=8, alpha=0.7)
                
                ax4.set_xlabel('TPSA')
                ax4.set_ylabel('Viscosity (Pa·s)')
                ax4.set_title('TPSA vs Viscosity')
                ax4.grid(True, alpha=0.3)
        
        # 调整布局
        plt.subplots_adjust(left=0.1, right=0.9, top=0.92, bottom=0.1, 
                          wspace=0.3, hspace=0.4)
        
        plt.savefig('Figure_Molecular_Correlations.png', dpi=300, bbox_inches='tight')
        print("   ✓ Figure saved: Figure_Molecular_Correlations.png")
        
        # 显示关键相关性
        if 'RPI' in self.parameters_df.columns:
            print("\n   Key Correlations with RPI:")
            for mol_prop in available_molecular:
                clean_data = self.parameters_df[[mol_prop, 'RPI']].dropna()
                if len(clean_data) >= 3:
                    corr, p_val = stats.pearsonr(clean_data[mol_prop], clean_data['RPI'])
                    significance = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
                    print(f"     {mol_prop}: r = {corr:.3f}{significance} (p = {p_val:.4f})")
        
        # 执行回归分析
        self._perform_regression_analysis(available_molecular, available_rheology)
        
        return fig
    
    def _perform_regression_analysis(self, molecular_params, rheology_params):
        """执行回归分析"""
        print("\n   Regression Analysis Results:")
        print("   " + "-" * 50)
        
        if 'RPI' not in rheology_params:
            return
        
        # 准备数据
        X_data = []
        y_data = []
        valid_salts = []
        
        for idx, row in self.parameters_df.iterrows():
            salt = row['盐类']
            # 检查是否有所有需要的分子性质
            has_all_params = True
            x_values = []
            
            for param in molecular_params:
                if param in row.index and not pd.isna(row[param]):
                    x_values.append(row[param])
                else:
                    has_all_params = False
                    break
            
            if has_all_params and not pd.isna(row['RPI']):
                X_data.append(x_values)
                y_data.append(row['RPI'])
                valid_salts.append(salt)
        
        if len(X_data) < 3:
            print("   Insufficient data for regression analysis")
            return
        
        X = np.array(X_data)
        y = np.array(y_data)
        
        # 执行多元线性回归
        model = LinearRegression()
        model.fit(X, y)
        
        print(f"   Samples used: {len(valid_salts)} salts")
        print(f"   R² score: {model.score(X, y):.4f}")
        
        # 打印系数
        print("\n   Regression Coefficients:")
        for i, param in enumerate(molecular_params):
            print(f"     {param}: {model.coef_[i]:.4f}")
        print(f"     Intercept: {model.intercept_:.4f}")
    
    def save_all_results(self):
        """保存所有结果，并导出绘图原始数据"""
        print("\n" + "=" * 80)
        print("[STEP 8] Saving all results...")
        
        os.makedirs('NatureComms_Results', exist_ok=True)
        
        # ---------- 原有的保存功能 ----------
        # 保存参数数据
        if not self.parameters_df.empty:
            csv_path = 'NatureComms_Results/Supplementary_Table_S1_Parameters.csv'
            self.parameters_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f"   ✓ Supplementary Table S1: {csv_path}")
            
            try:
                excel_path = 'NatureComms_Results/Supplementary_Table_S1_Parameters.xlsx'
                self.parameters_df.to_excel(excel_path, index=False)
                print(f"   ✓ Supplementary Table S1 (Excel): {excel_path}")
            except Exception as e:
                print(f"   ⚠ Excel save failed: {e}")
        
        # 保存分子性质数据
        if self.molecular_properties:
            mol_df = pd.DataFrame.from_dict(self.molecular_properties, orient='index')
            mol_df.index.name = 'Salt'
            mol_df.to_csv('NatureComms_Results/Molecular_Properties.csv')
            print(f"   ✓ Molecular properties: NatureComms_Results/Molecular_Properties.csv")
        
        # 移动和重命名图形
        figure_mapping = {
            'Figure_Salt_Comparison.png': 'NatureComms_Results/Fig1_Salt_Comparison.png',
            'Figure_2D_PCA_Analysis.png': 'NatureComms_Results/Fig2_PCA_Analysis.png',
            'Figure_Feature_Cluster_Analysis.png': 'NatureComms_Results/Fig3_Feature_Cluster.png',
            'Figure_3D_Parameter_Space.png': 'NatureComms_Results/Fig4_3D_Space.png',
            'Figure_Parallel_Coordinates.png': 'NatureComms_Results/Fig5_Parallel_Coordinates.png',
            'Figure_Molecular_Correlations.png': 'NatureComms_Results/Fig6_Molecular_Correlations.png',
        }
        
        for old_name, new_name in figure_mapping.items():
            if os.path.exists(old_name):
                try:
                    os.rename(old_name, new_name)
                    print(f"   ✓ {new_name}")
                except Exception as e:
                    print(f"   ⚠ Error renaming {old_name}: {e}")
        
        # ---------- 新增：导出用于外部绘图的原始数据 ----------
        print("\n   Exporting raw data for external plotting...")
        
        if not self.parameters_df.empty:
            # 1. 完整数据集（含所有参数、RPI、PCA坐标、聚类标签、分子性质）
            full_data = self.parameters_df.copy()
            # 确保 PCA 坐标和聚类标签存在（在 create_2d_network_analysis 中已添加）
            if 'PC1' not in full_data.columns:
                print("   Warning: PC1/PC2 not found, PCA may not have been run.")
            full_data_path = 'NatureComms_Results/Full_Data_for_Plotting.csv'
            full_data.to_csv(full_data_path, index=False, encoding='utf-8-sig')
            print(f"   ✓ Full dataset (including PCA coordinates): {full_data_path}")
            
            # 2. 专门导出 PCA 坐标（简化版）
            if 'PC1' in full_data.columns and 'PC2' in full_data.columns:
                pca_out = full_data[['盐类', 'RPI', 'PC1', 'PC2']].copy()
                if 'Cluster' in full_data.columns:
                    pca_out['Cluster'] = full_data['Cluster']
                pca_path = 'NatureComms_Results/PCA_Coordinates.csv'
                pca_out.to_csv(pca_path, index=False, encoding='utf-8-sig')
                print(f"   ✓ PCA coordinates: {pca_path}")
            
            # 3. 分子性质与流变参数的合并数据
            mol_rheo_cols = []
            # 分子性质列（常见）
            mol_cols = ['LogP', 'TPSA', 'MW', 'HBD', 'HBA', 'RotBonds', 'Charge', 'IonicRadius', 'HydrationEnergy']
            # 流变关键列
            rheo_cols = ['G_prime_linear', 'avg_viscosity', 'tan_delta_avg', 'Payne_effect', 
                         'recovery_percentage', 'Strain_stiffening', 'flow_index', 'consistency_index', 'RPI']
            
            selected_cols = ['盐类']
            for col in mol_cols + rheo_cols:
                if col in full_data.columns:
                    selected_cols.append(col)
            
            if len(selected_cols) > 1:
                mol_rheo_df = full_data[selected_cols].copy()
                mol_rheo_path = 'NatureComms_Results/Molecular_Rheology_Data.csv'
                mol_rheo_df.to_csv(mol_rheo_path, index=False, encoding='utf-8-sig')
                print(f"   ✓ Molecular-rheology combined data: {mol_rheo_path}")
        
        # 保存分析报告（原有）
        self._save_analysis_report()
        
        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        
        # 显示关键发现
        if not self.parameters_df.empty and 'RPI' in self.parameters_df.columns:
            print("\nKEY FINDINGS:")
            print("-" * 40)
            
            top5 = self.parameters_df.sort_values('RPI', ascending=False).head(5)
            print(f"\nTop 5 Performing Salts:")
            for i, (_, row) in enumerate(top5.iterrows(), 1):
                print(f"  {i}. {row['盐类']} (RPI: {row['RPI']:.3f})")
            
            if 'LogP' in self.parameters_df.columns:
                clean_data = self.parameters_df[['LogP', 'RPI']].dropna()
                if len(clean_data) >= 3:
                    corr_logp, p_logp = stats.pearsonr(clean_data['LogP'], clean_data['RPI'])
                    print(f"\nLogP vs RPI correlation: r = {corr_logp:.3f} (p = {p_logp:.4f})")
            
            if 'TPSA' in self.parameters_df.columns:
                clean_data = self.parameters_df[['TPSA', 'RPI']].dropna()
                if len(clean_data) >= 3:
                    corr_tpsa, p_tpsa = stats.pearsonr(clean_data['TPSA'], clean_data['RPI'])
                    print(f"TPSA vs RPI correlation: r = {corr_tpsa:.3f} (p = {p_tpsa:.4f})")
    
    def _save_analysis_report(self):
        """保存详细分析报告"""
        if self.parameters_df.empty:
            return
        
        report_path = 'NatureComms_Results/Analysis_Report.txt'
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("RHEOLOGICAL ANALYSIS REPORT\n")
            f.write("Nature Communications Submission\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Report Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("1. STUDY OVERVIEW\n")
            f.write("-" * 40 + "\n")
            f.write(f"Number of salts analyzed: {len(self.parameters_df)}\n")
            f.write(f"Number of parameters extracted: {len(self.parameters_df.columns)}\n")
            f.write(f"Molecular properties included: {len(self.molecular_properties)} salts\n\n")
            
            f.write("2. KEY RESULTS\n")
            f.write("-" * 40 + "\n")
            
            if 'RPI' in self.parameters_df.columns:
                f.write("2.1 Rheological Performance Ranking:\n")
                ranked = self.parameters_df.sort_values('RPI', ascending=False)
                for i, (_, row) in enumerate(ranked.iterrows(), 1):
                    f.write(f"  {i:2d}. {row['盐类']:10s} RPI = {row['RPI']:.3f}\n")
                f.write("\n")
            
            f.write("3. MOLECULAR PROPERTIES ANALYSIS\n")
            f.write("-" * 40 + "\n")
            
            if self.molecular_properties:
                f.write("3.1 Molecular Properties Loaded:\n")
                for salt, props in self.molecular_properties.items():
                    f.write(f"  {salt}: {', '.join(props.keys())}\n")
                f.write("\n")
            
            f.write("4. METHODOLOGY\n")
            f.write("-" * 40 + "\n")
            f.write("4.1 Data Processing:\n")
            f.write("  - Multi-sheet Excel/CSV data loading\n")
            f.write("  - Rheological parameter extraction\n")
            f.write("  - Molecular property integration\n\n")
            
            f.write("4.2 Analysis Methods:\n")
            f.write("  - Principal Component Analysis (PCA)\n")
            f.write("  - K-means clustering\n")
            f.write("  - Correlation analysis\n")
            f.write("  - 2D/3D visualization\n\n")
            
            f.write("5. DATA AVAILABILITY\n")
            f.write("-" * 40 + "\n")
            f.write("All data are available in the supplementary files.\n\n")
            
            f.write("=" * 80 + "\n")
            f.write("END OF REPORT\n")
            f.write("=" * 80 + "\n")
        
        print(f"   ✓ Analysis report: {report_path}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================
def main():
    """主函数"""
    try:
        analyzer = NatureCommsRheologyAnalyzer(data_path='F:/Data/rheology_data')
        
        print("\n" + "=" * 80)
        print("NATURE COMMUNICATIONS RHEOLOGICAL ANALYSIS")
        print("=" * 80)
        
        # 1. 加载所有sheet数据
        if not analyzer.load_all_sheets():
            print("Data loading failed, please check file path and format")
            return
        
        # 2. 按盐类合并数据
        analyzer.combine_data_by_salt()
        
        # 3. 提取流变学参数
        analyzer.extract_rheological_parameters()
        
        # 4. 分析盐类效应
        analyzer.analyze_salt_effects()
        
        # 5. 创建二维网络分析
        analyzer.create_2d_network_analysis()
        
        # 6. 创建三维可视化
        analyzer.create_3d_visualization()
        
        # 7. 创建分子性质相关性分析
        analyzer.create_molecular_correlation_analysis()
        
        # 8. 保存所有结果（包括导出的原始数据）
        analyzer.save_all_results()
        
    except Exception as e:
        print(f"\nERROR: {e}")
        print("\nTraceback:")
        print(traceback.format_exc())

if __name__ == "__main__":
    main()