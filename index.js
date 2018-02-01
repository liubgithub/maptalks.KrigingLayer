import * as maptalks from 'maptalks';
import K from 'kriging';

/*
width is the pixel cell's width
model is the kriging render method,mainly includ 'Gaussian','Exponential',and 'Spherical'
Gaussian: k(a,b) = w[0] + w[1] * ( 1 - exp{ -( ||a-b|| / range )2 / A } )
Exponential: k(a,b) = w[0] + w[1] * ( 1 - exp{ -( ||a-b|| / range ) / A } )
Spherical: k(a,b) = w[0] + w[1] * ( 1.5 * ( ||a-b|| / range ) - 0.5 * ( ||a-b|| / range )3 )
*/
const options = {
    width:0.001,
    model:'exponential',
    sigma2:0,
    alpha:10
};

export class KrigingLayer extends maptalks.Layer {

    constructor(id, interest, options) {
        if (!Array.isArray(interest)) {
            options = interest;
            interest = null;
        }
        super(id, options);
        this._interest = interest || [];
    }

    getData() {
        return this._interest;
    }

    setData(interest) {
        this._interest = interest || [];
        return this.redraw();
    }

    onConfig(conf) {
        super.onConfig.apply(this, arguments);
        for (const p in conf) {
            if (options[p]) {
                return this.redraw();
            }
        }
        return this;
    }

    setModel(model) {
        if (model instanceof String) {
            this.options['model'] = model;
            return this.redraw();
        }
        return this;
    }

    redraw() {
        const renderer = this._getRenderer();
        if (renderer) {
            renderer.clearHeatCache();
            renderer.setToRedraw();
        }
        return this;
    }

    isEmpty() {
        if (!this._interest.length) {
            return true;
        }
        return false;
    }

    clear() {
        this._interest = [];
        this.redraw();
        this.fire('clear');
        return this;
    }

    /**
     * Export the KrigingLayer's JSON.
     * @return {Object} layer's JSON
     */
    toJSON(options) {
        if (!options) {
            options = {};
        }
        const json = {
            'type'      : this.getJSONType(),
            'id'        : this.getId(),
            'options'   : this.config()
        };
        const data = this.getData();
        if (options['clipExtent']) {
            let clipExtent = new maptalks.Extent(options['clipExtent']);
            const r = this._getHeatRadius();
            if (r) {
                clipExtent = clipExtent._expand(r);
            }
            const clipped = [];
            for (let i = 0, len = data.length; i < len; i++) {
                if (clipExtent.contains(new maptalks.Coordinate(data[i][0], data[i][1]))) {
                    clipped.push(data[i]);
                }
            }
            json['data'] = clipped;
        } else {
            json['data'] = data;
        }

        return json;
    }

    /**
     * Reproduce a KrigingLayer from layer's JSON.
     * @param  {Object} json - layer's JSON
     * @return {maptalks.KrigingLayer}
     * @static
     * @private
     * @function
     */
    static fromJSON(json) {
        if (!json || json['type'] !== 'KrigingLayer') {
            return null;
        }
        return new KrigingLayer(json['id'], json['data'], json['options']);
    }
}

KrigingLayer.mergeOptions(options);

KrigingLayer.registerJSONType('KrigingLayer');

KrigingLayer.registerRenderer('canvas', class extends maptalks.renderer.CanvasRenderer {

    needToRedraw() {
        const map = this.layer.getMap();
        if (map.isZooming()) {
            return false;
        }
        if (map.isMoving()) {
            return false;
        }
        return super.needToRedraw();
    }

    draw() {
        this.prepareCanvas();
        if (!this._isInExtent()) {
            this.completeRender();
            return;
        }
        this._plot();
        this.completeRender();
    }

    drawOnInteracting() {
        this.draw();
    }

    _plot() {
        const map = this.layer.getMap();
        const width = this.layer.options['width'];
        const colors = this.layer.options['colors'];
        const regions = this.layer.options['regions'];
        const model = this.layer.options['model'];
        const sigma2 = this.layer.options['sigma2'];
        const alpha = this.layer.options['alpha'];
        const _polygons = this._handRegions(regions);
        const extent = map.getExtent();
        const data = this.layer.getData();
        const lngs = data.map(function (d) {
            return d[0];
        });
        const lats = data.map(function (d) {
            return d[1];
        });
        const values = data.map(function (d) {
            return d[2];
        });
        const variogram = K.kriging.train(values, lngs, lats, model, sigma2, alpha);
        const grid = K.kriging.grid(_polygons, variogram, width);
        K.kriging.plot(this.canvas, grid, [extent.xmin, extent.xmax], [extent.ymin, extent.ymax], colors);
    }

    _isInExtent() {
        const map = this.layer.getMap();
        const mapExtent = map.getExtent();
        const regions = this.layer.options['regions'];
        const regionExtent = regions.getExtent();
        if (mapExtent.intersects(regionExtent))
            return true;
        else return false;
    }

    onZoomEnd() {
        //delete this._heatViews;
        super.onZoomEnd.apply(this, arguments);
    }

    onResize() {
        //this._interest._width  = this.canvas.width;
        //this._interest._height = this.canvas.height;
        super.onResize.apply(this, arguments);
    }

    onDragRotateEnd(e) {
        super.onDragRotateEnd(e);
    }

    resizeCanvas() {
        if (!this.canvas) {
            return;
        }
    }

    clearCanvas() {
        if (!this.canvas) {
            return;
        }

    }

    _handRegions(regions) {
        const _polygons = regions.getCoordinates().map(function (coords) {
            return coords.map(function (c) {
                return [c.x, c.y];
            });
        });
        return _polygons;
    }
});
