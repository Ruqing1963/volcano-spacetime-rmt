"""Export the four target datasets to tidy CSV files for the repository."""
import csv
import numpy as np
import importlib.util

spec = importlib.util.spec_from_file_location("vp", "/home/claude/volcano_rmt_pipeline.py")
vp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vp)

OUT = "/home/claude/repo_build/volcano-spacetime-rmt/data"

# ---- Target A: Etna flank eruptions ----
etna = vp.load_etna_eruptions()
with open(f"{OUT}/target_A_etna_flank_eruptions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["index", "eruption_year_CE"])
    for i, y in enumerate(np.sort(etna), 1):
        w.writerow([i, int(y)])
print(f"A: {len(etna)} Etna eruptions")

# ---- Target B: six-volcano pooled source eruption catalogues ----
import inspect
src = inspect.getsource(vp.load_multi_volcano_pooled)
# Reconstruct the volcano dict by executing the loader and grabbing dates:
# Easiest: re-declare the same six arrays by calling a helper exposed below.
volcanoes = {
    'Mauna Loa (Hawaii hotspot)': [
        1843,1849,1851,1852,1855,1859,1868,1871,1873,1875,1876,1877,1879,1880,
        1887,1892,1896,1899,1903,1907,1914,1916,1919,1926,1933,1935,1940,1942,
        1949,1950,1975,1984,2022],
    'Piton de la Fournaise (Reunion hotspot)': [
        1640,1649,1670,1708,1717,1721,1733,1751,1753,1759,1766,1774,1775,1776,
        1784,1786,1787,1789,1791,1797,1800,1802,1812,1813,1814,1820,1824,1830,
        1832,1843,1844,1849,1852,1858,1859,1860,1861,1863,1868,1869,1871,1872,
        1874,1878,1882,1894,1898,1901,1903,1907,1909,1913,1915,1920,1924,1926,
        1929,1930,1931,1932,1933,1934,1935,1936,1937,1938,1941,1942,1943,1944,
        1945,1947,1948,1950,1951,1953,1954,1955,1956,1957,1958,1959,1960,1961,
        1963,1964,1966,1972,1973,1975,1976,1977,1979,1981,1983,1984,1985,1986,
        1987,1988,1990,1991,1992,1998,2000,2001,2002,2003,2004,2005,2006,2007,
        2008,2009,2010,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024],
    'Hekla (Iceland rift)': [
        1104,1158,1206,1222,1300,1341,1389,1510,1597,1636,1693,1725,1766,1845,
        1878,1913,1947,1970,1980,1981,1991,2000],
    'Vesuvius (subduction, Italy)': [
        79,203,379,472,512,536,685,787,860,968,991,999,1007,1037,1068,1078,
        1139,1150,1270,1347,1500,1631,1637,1649,1654,1660,1682,1694,1698,1707,
        1714,1717,1723,1730,1737,1751,1754,1760,1767,1770,1779,1785,1790,1794,
        1804,1806,1810,1813,1817,1822,1831,1833,1834,1839,1841,1845,1850,1855,
        1858,1861,1868,1871,1872,1875,1891,1895,1899,1906,1926,1929,1944],
    'Klyuchevskoy (Kamchatka subduction)': [
        1697,1708,1710,1711,1720,1725,1727,1729,1731,1737,1740,1746,1756,1757,
        1762,1767,1769,1772,1778,1785,1787,1788,1789,1790,1791,1793,1795,1797,
        1811,1813,1820,1822,1824,1825,1829,1831,1835,1840,1841,1848,1852,1853,
        1854,1855,1856,1857,1858,1859,1861,1862,1865,1868,1870,1873,1877,1878,
        1879,1882,1883,1885,1886,1887,1890,1892,1893,1894,1895,1896,1897,1898,
        1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,1911,1912,1913,
        1914,1915,1919,1922,1923,1925,1926,1928,1931,1932,1935,1937,1938,1939,
        1943,1944,1945,1946,1948,1951,1953,1954,1956,1957,1958,1960,1961,1962,
        1963,1965,1966,1967,1968,1969,1970,1971,1972,1974,1977,1978,1979,1980,
        1981,1982,1983,1984,1985,1986,1987,1988,1989,1990,1993,1994,1996,1997,
        1998,2002,2003,2004,2005,2007,2008,2009,2010,2012,2013,2015,2016,2017,
        2018,2019,2020,2021,2022,2023],
    'Colima (Mexico subduction)': [
        1560,1576,1585,1590,1606,1611,1613,1622,1680,1690,1711,1749,1770,1795,
        1806,1818,1869,1872,1878,1880,1885,1889,1890,1893,1903,1908,1913,1917,
        1925,1931,1941,1957,1961,1975,1981,1982,1985,1987,1991,1994,1998,1999,
        2001,2004,2005,2007,2014,2015,2016],
}
with open(f"{OUT}/target_B_multivolcano_eruptions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["volcano", "tectonic_setting", "eruption_year_CE"])
    for name, yrs in volcanoes.items():
        base = name.split(" (")[0]
        setting = name.split("(")[1].rstrip(")") if "(" in name else ""
        for y in sorted(yrs):
            w.writerow([base, setting, int(y)])
nB = sum(len(v) for v in volcanoes.values())
print(f"B: {nB} eruptions across {len(volcanoes)} volcanoes")

# ---- Target C: Cascade arc positions ----
casc = vp.load_cascade_arc()
with open(f"{OUT}/target_C_cascade_arc_positions.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["volcano", "along_arc_distance_km"])
    for name, d in sorted(casc.items(), key=lambda kv: kv[1]):
        w.writerow([name, d])
print(f"C: {len(casc)} Cascade volcanoes")

# ---- Target D: Iceland vents ----
ice = vp.load_iceland_rift()
positions = ice["positions"]
# Reconstruct cluster membership for documentation
clusters = {
    "Reykjanes": [0,1.5,3,5,6,8,11,14,16,19,22,26],
    "Hengill": [48,49,50.5,52,53,55,57,60],
    "Eyjafjallajokull-Katla": [76,78,79,81,83,85,87,89,92,95],
    "Hekla-Torfajokull": [113,114.5,116,118,120,122],
    "Vatnajokull": [152,153.5,155,156,157.5,159,160,161,163,165,167,169,171,173,175],
    "Askja-Herdubreid": [198,200,201,203,205,207,210],
    "Krafla-Theistareykir": [238,240,241,243,245,248,250,254,258],
}
pos2cluster = {}
for cl, vals in clusters.items():
    for v in vals:
        pos2cluster[round(v, 2)] = cl
with open(f"{OUT}/target_D_iceland_vents.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["index", "transect_position_km", "volcanic_system"])
    for i, p in enumerate(np.sort(positions), 1):
        w.writerow([i, round(float(p), 2), pos2cluster.get(round(float(p), 2), "")])
print(f"D: {len(positions)} Iceland vents")
print("All CSVs written to", OUT)
