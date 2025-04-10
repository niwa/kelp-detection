# reference from https://www.dash-leaflet.com/

from dash import Dash
import dash_leaflet as dl

app = Dash()
'''app.layout = dl.Map(dl.TileLayer(), center=[56,10], zoom=6, style={'height': '50vh'})'''

'''app.layout = dl.Map([
    dl.TileLayer(),
    dl.WMSTileLayer(url="https://mesonet.agron.iastate.edu/cgi-bin/wms/nexrad/n0r.cgi",
                    layers="nexrad-n0r-900913", format="image/png", transparent=True)
], center=[40, -100], zoom=4, style={'height': '50vh'})'''

app.layout = dl.Map([
    dl.TileLayer(),
    dl.WMSTileLayer(url="https://data.linz.govt.nz/services;key=62319684538a43d7b3a63cbebf678ee4/wmts/1.0.0/layer/109401/WMTSCapabilities.xml",
                    layers="109401") # 109401-nz-10m-satellite-imagery-2021-2022
], center=[40, -100], zoom=4, style={'height': '50vh'})
# https://basemaps.linz.govt.nz/v1/tiles/aerial/WebMercatorQuad/WMTSCapabilities.xml?api=c01jjaz37wdwavavyw3q2gtrbdp

if __name__ == '__main__':
    app.run_server()