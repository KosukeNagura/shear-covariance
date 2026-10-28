import numpy as np
from scipy.special import eval_legendre
import healpy as hp

def calculate_2pcf_covariance(pseudo_cl_cov, n_a, n_b, n_c, n_d, theta_bins, N_pair):
    """
    論文の式(67)を用いて、2点相関関数の共分散行列を計算する。
    
    Parameters:
    -----------
    pseudo_cl_cov : numpy.ndarray
        pseudo-Clの共分散行列。形状は (ell_max, ell_max)。
    n_a, n_b, n_c, n_d : float
        各トレーサーの数密度。
    theta_bins : numpy.ndarray
        thetaのビン境界（ラジアン）。長さは num_bins + 1。
    N_pair : numpy.ndarray
        各thetaビンにおける観測されたペア数。長さは num_bins。
        
    Returns:
    --------
    xi_cov : numpy.ndarray
        2点相関関数の共分散行列。形状は (num_bins, num_bins)。
    """
    ell_max = pseudo_cl_cov.shape[0] - 1
    num_bins = len(theta_bins) - 1
    
    # thetaの境界をx = cos(theta)に変換
    # 式(67)の積分範囲に対応 [1, 2]
    x_bins = np.cos(theta_bins)
    
    # Delta P_ell行列の作成: 形状 (num_bins, ell_max + 1)
    # delta_p[i, ell] = [P_{ell+1}(x) - P_{ell-1}(x)]_{theta_i-}^{theta_i+}
    delta_p = np.zeros((num_bins, ell_max + 1))
    
    for i in range(num_bins):
        # 積分範囲のx座標（cos(theta)なので、theta_minがx_highになる）
        x_high = x_bins[i]   # theta_ab-
        x_low = x_bins[i+1]  # theta_ab+
        
        for ell in range(ell_max + 1):
            # P_{ell+1}の差分
            term_plus = eval_legendre(ell + 1, x_high) - eval_legendre(ell + 1, x_low)
            
            # P_{ell-1}の差分（ell=0の場合は0とする）
            if ell > 0:
                term_minus = eval_legendre(ell - 1, x_high) - eval_legendre(ell - 1, x_low)
            else:
                # 論文の記述に基づき、負の次数は寄与しない [2]
                term_minus = 0
            
            # 修正前
            delta_p[i, ell] = term_plus - term_minus

            # 修正後: 必ず (2*ell + 1) で割る → やっぱこっちが間違い
            # if ell > 0:
                # delta_p[i, ell] = (term_plus - term_minus) / (2 * ell + 1)
            # else:
            # ell=0 の積分は単に区間長 (x_high - x_low)
            # P_1(x) - P_{-1}(x) = x - (-x)? 定義によるが、通常 ell=0 は別処理したほうが安全
            # 現在のコードフローなら以下でOK (term_minusがどうなるかによるが、明示的に書くのが安全)
                # delta_p[i, ell] = x_high - x_low

    # (2 * pi)^2 の定数項 [1]
    const_factor = (2 * np.pi)**2
    
    # 行列演算による2重和 \sum_{l1, l2} の実行
    # sum_part[i, j] = \sum_{l1, l2} DeltaP[i, l1] * Cov(l1, l2) * DeltaP[j, l2]
    sum_part = delta_p @ pseudo_cl_cov @ delta_p.T
    
    # N_pairによる正規化項の作成
    # 1 / (N_pair_i * N_pair_j) の外積行列
    inv_n_pair = 1.0 / N_pair
    norm_matrix = np.outer(inv_n_pair, inv_n_pair)
    
    # 全体の数密度因子 [1]
    density_factor = n_a * n_b * n_c * n_d
    
    # 式(67)を解いてxi_covを算出
    # Cov(xi) = (n_a*n_b*n_c*n_d / (N_p_i * N_p_j)) * (2pi)^2 * sum_part
    xi_cov = density_factor * norm_matrix * const_factor * sum_part
    
    return xi_cov





def calculate_number_density(count_map, weight_map):
    """
    HEALPixマップから、重みを考慮した銀河の平均数密度 n を計算する。
    
    Parameters:
    -----------
    count_map : numpy.ndarray
        ピクセルごとの銀河数が入ったHEALPixマップ。
    weight_map : numpy.ndarray
        各ピクセルの重み（完全性やマスク）が入ったマップ。
        未観測領域は0、完全にカバーされている領域は1など。
        
    Returns:
    --------
    n_bar : float
        単位立体角（steradian）あたりの平均数密度。
    """
    # マップの解像度（nside）を取得
    nside = hp.get_nside(count_map)
    
    # 1ピクセルあたりの面積（steradians）を取得 [2]
    omega_pix = hp.nside2pixarea(nside)
    
    # マスクが適用されている（weight > 0）ピクセルのみを対象にする
    mask = weight_map > 0
    
    # ピクセルあたりの平均ソース数 N_bar を計算 [1, 2]
    # N_bar = sum(N_i) / sum(W_i)
    total_counts = np.sum(count_map[mask])
    total_weights = np.sum(weight_map[mask])
    
    n_per_pixel = total_counts / total_weights
    
    # 立体角あたりの数密度 n に変換 [2]
    # n_bar = N_bar / omega_pix
    n_bar = n_per_pixel / omega_pix
    
    return n_bar




def calculate_number_density_r(count_map, weight_map):
    """
    HEALPixマップから、重みを考慮した銀河の平均数密度 n を計算する。
    
    Parameters:
    -----------
    count_map : numpy.ndarray
        ピクセルごとの銀河数が入ったHEALPixマップ。
    weight_map : numpy.ndarray
        各ピクセルの重み（完全性やマスク）が入ったマップ。
        未観測領域は0、完全にカバーされている領域は1など。
        
    Returns:
    --------
    n_bar : float
        単位立体角（steradian）あたりの平均数密度。
    """
    # マップの解像度（nside）を取得
    nside = hp.get_nside(count_map)
    
    # 1ピクセルあたりの面積（steradians）を取得 [2]
    omega_pix = hp.nside2pixarea(nside)
    
    # ピクセルあたりの平均ソース数 N_bar を計算 [1, 2]
    # N_bar = sum(N_i) / sum(W_i)
    total_counts = np.sum(count_map)
    total_weights = np.sum(weight_map)
    
    n_per_pixel = total_counts / total_weights
    
    # 立体角あたりの数密度 n に変換 [2]
    # n_bar = N_bar / omega_pix
    n_bar = n_per_pixel / omega_pix
    
    return n_bar





def galaxy_catalog_to_healpix(
    ra_deg,
    dec_deg,
    nside,
    weight=None,
    nest=False,
    return_delta=False
):
    """
    Parameters
    ----------
    ra_deg : array_like
        Right Ascension in degrees [0, 360)
    dec_deg : array_like
        Declination in degrees [-90, 90]
    nside : int
        HEALPix nside
    weight : array_like or None
        Optional galaxy weight
    nest : bool
        Use NESTED ordering if True, else RING
    return_delta : bool
        If True, also return overdensity map

    Returns
    -------
    count_map : ndarray
        Galaxy number count per pixel
    weighted_map : ndarray or None
        Weighted map (None if weight is None)
    delta_map : ndarray or None
        Overdensity map (None if return_delta=False)
    """

    ra = np.asarray(ra_deg)
    dec = np.asarray(dec_deg)

    # healpy uses theta = colatitude
    theta = np.radians(90.0 - dec)
    phi   = np.radians(ra)

    npix = hp.nside2npix(nside)

    # pixel indices
    ipix = hp.ang2pix(nside, theta, phi, nest=nest)

    # number count map
    count_map = np.zeros(npix, dtype=float)
    np.add.at(count_map, ipix, 1.0)

    weighted_map = None
    if weight is not None:
        weight = np.asarray(weight)
        weighted_map = np.zeros(npix, dtype=float)
        np.add.at(weighted_map, ipix, weight)

    delta_map = None
    if return_delta:
        mean = count_map.mean()
        delta_map = np.zeros_like(count_map)
        if mean > 0:
            delta_map = count_map / mean - 1.0

    return count_map, weighted_map, delta_map


def calculate_n_pair(weight_map, n_a, n_b, theta_bins):
    """
    HEALPixマスクマップから、各thetaビンにおける期待ペア数 N_pair を計算する。
    論文の式(C5)および付録Cの導出に基づく。
    
    Parameters:
    -----------
    weight_map : numpy.ndarray
        サーベイのマスク（完全性）が入ったHEALPixマップ。
    n_a, n_b : float
        各トレーサーの単位立体角（steradian）あたりの数密度。
    theta_bins : numpy.ndarray
        thetaのビン境界（ラジアン）。長さは num_bins + 1。
        
    Returns:
    --------
    n_pair : numpy.ndarray
        各角度ビンにおける期待ペア数（1D配列）。
    """
    # 1. マスクの角度パワースペクトル C_W_ell を計算
    # hp.anafast は規格化された C_ell (unit: steradian相当) を返す
    cl_w = hp.anafast(weight_map)
    ell_max = len(cl_w) - 1
    num_bins = len(theta_bins) - 1
    
    # 2. thetaの境界を x = cos(theta) に変換
    # 積分範囲 [theta_min, theta_max] は [x_high, x_low] に対応
    x_bins = np.cos(theta_bins)
    
    n_pair_list = []
    
    # 各ビンについて積分を実行
    for i in range(num_bins):
        x_high = x_bins[i]
        x_low = x_bins[i+1]
        
        # 3. 式(C5)の総和部分の計算
        # ell=0 の項: Integral(P_0)dx = x
        #sum_val = cl_w * (x_high - x_low)  こっちミス
        sum_val = cl_w[0] * (x_high - x_low)
        
        # ell >= 1 の項: Integral(P_ell)dx = [P_{ell+1} - P_{ell-1}] / (2*ell + 1)
        # ただし式(C5)の表記に合わせ、分母の(2*ell+1)が相殺された形を実装
        ells = np.arange(1, ell_max + 1)
        
        # P_{ell+1}(x) - P_{ell-1}(x) の境界値の差
        delta_p_high = eval_legendre(ells + 1, x_high) - eval_legendre(ells - 1, x_high)
        delta_p_low = eval_legendre(ells + 1, x_low) - eval_legendre(ells - 1, x_low)
        delta_p = delta_p_high - delta_p_low
        
        # C_W_ell * Delta_P の和をとる
        sum_val += np.sum(cl_w[1:] * delta_p)

        
        # 4. 全体の定数倍を適用
        # 式(C5)より: (8 * pi^2 * n_a * n_b) * (1 / 4 * pi) * sum_val
        # = 2 * pi * n_a * n_b * sum_val
        n_pair_val = 2 * np.pi * n_a * n_b * sum_val
        n_pair_list.append(n_pair_val)
        
    return np.array(n_pair_list)

# n_pair_expected = calculate_n_pair(weight_map, n_dens, n_dens, theta_edges)