import * as maptalks from 'maptalks';
import K from 'kriging';

const options = {
    'max' : 1,
    'gradient' : {
        0.4: 'blue',
        0.6: 'cyan',
        0.7: 'lime',
        0.8: 'yellow',
        1.0: 'red'
    },
    'radius' : 25,
    'blur' : 15,
    'minOpacity' : 0.05
};

export class KrigingLayer extends maptalks.Layer {

    constructor(id, heats, options) {
        if (!Array.isArray(heats)) {
            options = heats;
            heats = null;
        }
        super(id, options);
        this._heats = heats || [];
    }

    getData() {
        return this._heats;
    }

    setData(heats) {
        this._heats = heats || [];
        return this.redraw();
    }

    addPoint(heat) {
        if (!heat) {
            return this;
        }
        if (heat[0] && Array.isArray(heat[0])) {
            maptalks.Util.pushIn(this._heats, heat);
        } else {
            this._heats.push(heat);
        }
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

    redraw() {
        const renderer = this._getRenderer();
        if (renderer) {
            renderer.clearHeatCache();
            renderer.setToRedraw();
        }
        return this;
    }

    isEmpty() {
        if (!this._heats.length) {
            return true;
        }
        return false;
    }

    clear() {
        this._heats = [];
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
        if (!json || json['type'] !== 'KrigingLayer') { return null; }
        return new KrigingLayer(json['id'], json['data'], json['options']);
    }


    _getHeatRadius() {
        if (!this._getRenderer()) {
            return null;
        }
        return this._getRenderer()._heatRadius;
    }
}

KrigingLayer.mergeOptions(options);

KrigingLayer.registerJSONType('KrigingLayer');

KrigingLayer.registerRenderer('canvas', class extends maptalks.renderer.CanvasRenderer {

    draw() {
        K.addd(1, 2);
        this.completeRender();
    }

    drawOnInteracting() {
        this.draw();
    }


    onZoomEnd() {
        delete this._heatViews;
        super.onZoomEnd.apply(this, arguments);
    }

    onResize() {
        this._heater._width  = this.canvas.width;
        this._heater._height = this.canvas.height;
        super.onResize.apply(this, arguments);
    }

    onRemove() {
        this.clearHeatCache();
        delete this._heater;
    }

    clearHeatCache() {
        delete this._heatViews;
    }
});
