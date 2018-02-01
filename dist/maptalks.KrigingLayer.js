/*!
 * maptalks.kriginglayer v0.1.0
 * LICENSE : MIT
 * (c) 2016-2018 maptalks.org
 */
/*!
 * requires maptalks@^0.37.0 
 */
(function (global, factory) {
	typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports, require('maptalks')) :
	typeof define === 'function' && define.amd ? define(['exports', 'maptalks'], factory) :
	(factory((global.maptalks = global.maptalks || {}),global.maptalks));
}(this, (function (exports,maptalks) { 'use strict';

// Extend the Array class
Array.prototype.max = function () {
	return Math.max.apply(null, this);
};
Array.prototype.min = function () {
	return Math.min.apply(null, this);
};
Array.prototype.mean = function () {
	var i = void 0,
	    sum = void 0;
	for (i = 0, sum = 0; i < this.length; i++) {
		sum += this[i];
	}return sum / this.length;
};
Array.prototype.rep = function (n) {
	return Array.apply(null, new Array(n)).map(Number.prototype.valueOf, this[0]);
};
Array.prototype.pip = function (x, y) {
	var i = void 0,
	    j = void 0,
	    c = false;
	for (i = 0, j = this.length - 1; i < this.length; j = i++) {
		if (this[i][1] > y != this[j][1] > y && x < (this[j][0] - this[i][0]) * (y - this[i][1]) / (this[j][1] - this[i][1]) + this[i][0]) {
			c = !c;
		}
	}
	return c;
};

var kriging = {};

// Matrix algebra
var kriging_matrix_diag = function kriging_matrix_diag(c, n) {
	var i = void 0,
	    Z = [0].rep(n * n);
	for (i = 0; i < n; i++) {
		Z[i * n + i] = c;
	}return Z;
};
var kriging_matrix_transpose = function kriging_matrix_transpose(X, n, m) {
	var i = void 0,
	    j = void 0,
	    Z = Array(m * n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			Z[j * n + i] = X[i * m + j];
		}
	}return Z;
};
var kriging_matrix_add = function kriging_matrix_add(X, Y, n, m) {
	var i = void 0,
	    j = void 0,
	    Z = Array(n * m);
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			Z[i * m + j] = X[i * m + j] + Y[i * m + j];
		}
	}return Z;
};
// Naive matrix multiplication
var kriging_matrix_multiply = function kriging_matrix_multiply(X, Y, n, m, p) {
	var i = void 0,
	    j = void 0,
	    k = void 0,
	    Z = Array(n * p);
	for (i = 0; i < n; i++) {
		for (j = 0; j < p; j++) {
			Z[i * p + j] = 0;
			for (k = 0; k < m; k++) {
				Z[i * p + j] += X[i * m + k] * Y[k * p + j];
			}
		}
	}
	return Z;
};
// Cholesky decomposition
var kriging_matrix_chol = function kriging_matrix_chol(X, n) {
	var i = void 0,
	    j = void 0,
	    k = void 0,
	    sum = void 0,
	    p = Array(n);
	for (i = 0; i < n; i++) {
		p[i] = X[i * n + i];
	}for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			p[i] -= X[i * n + j] * X[i * n + j];
		}if (p[i] <= 0) return false;
		p[i] = Math.sqrt(p[i]);
		for (j = i + 1; j < n; j++) {
			for (k = 0; k < i; k++) {
				X[j * n + i] -= X[j * n + k] * X[i * n + k];
			}X[j * n + i] /= p[i];
		}
	}
	for (i = 0; i < n; i++) {
		X[i * n + i] = p[i];
	}return true;
};
// Inversion of cholesky decomposition
var kriging_matrix_chol2inv = function kriging_matrix_chol2inv(X, n) {
	var i = void 0,
	    j = void 0,
	    k = void 0,
	    sum = void 0;
	for (i = 0; i < n; i++) {
		X[i * n + i] = 1 / X[i * n + i];
		for (j = i + 1; j < n; j++) {
			sum = 0;
			for (k = i; k < j; k++) {
				sum -= X[j * n + k] * X[k * n + i];
			}X[j * n + i] = sum / X[j * n + j];
		}
	}
	for (i = 0; i < n; i++) {
		for (j = i + 1; j < n; j++) {
			X[i * n + j] = 0;
		}
	}for (i = 0; i < n; i++) {
		X[i * n + i] *= X[i * n + i];
		for (k = i + 1; k < n; k++) {
			X[i * n + i] += X[k * n + i] * X[k * n + i];
		}for (j = i + 1; j < n; j++) {
			for (k = j; k < n; k++) {
				X[i * n + j] += X[k * n + i] * X[k * n + j];
			}
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			X[i * n + j] = X[j * n + i];
		}
	}
};
// Inversion via gauss-jordan elimination
var kriging_matrix_solve = function kriging_matrix_solve(X, n) {
	var m = n;
	var b = Array(n * n);
	var indxc = Array(n);
	var indxr = Array(n);
	var ipiv = Array(n);
	var i = void 0,
	    icol = void 0,
	    irow = void 0,
	    j = void 0,
	    k = void 0,
	    l = void 0,
	    ll = void 0;
	var big = void 0,
	    dum = void 0,
	    pivinv = void 0,
	    temp = void 0;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j) b[i * n + j] = 1;else b[i * n + j] = 0;
		}
	}for (j = 0; j < n; j++) {
		ipiv[j] = 0;
	}for (i = 0; i < n; i++) {
		big = 0;
		for (j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (Math.abs(X[j * n + k]) >= big) {
							big = Math.abs(X[j * n + k]);
							irow = j;
							icol = k;
						}
					}
				}
			}
		}
		++ipiv[icol];

		if (irow != icol) {
			for (l = 0; l < n; l++) {
				temp = X[irow * n + l];
				X[irow * n + l] = X[icol * n + l];
				X[icol * n + l] = temp;
			}
			for (l = 0; l < m; l++) {
				temp = b[irow * n + l];
				b[irow * n + l] = b[icol * n + l];
				b[icol * n + l] = temp;
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;

		if (X[icol * n + icol] == 0) return false; // Singular

		pivinv = 1 / X[icol * n + icol];
		X[icol * n + icol] = 1;
		for (l = 0; l < n; l++) {
			X[icol * n + l] *= pivinv;
		}for (l = 0; l < m; l++) {
			b[icol * n + l] *= pivinv;
		}for (ll = 0; ll < n; ll++) {
			if (ll != icol) {
				dum = X[ll * n + icol];
				X[ll * n + icol] = 0;
				for (l = 0; l < n; l++) {
					X[ll * n + l] -= X[icol * n + l] * dum;
				}for (l = 0; l < m; l++) {
					b[ll * n + l] -= b[icol * n + l] * dum;
				}
			}
		}
	}
	for (l = n - 1; l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for (k = 0; k < n; k++) {
				temp = X[k * n + indxr[l]];
				X[k * n + indxr[l]] = X[k * n + indxc[l]];
				X[k * n + indxc[l]] = temp;
			}
		}
	}return true;
};

// Variogram models
var kriging_variogram_gaussian = function kriging_variogram_gaussian(h, nugget, range, sill, A) {
	return nugget + (sill - nugget) / range * (1.0 - Math.exp(-(1.0 / A) * Math.pow(h / range, 2)));
};
var kriging_variogram_exponential = function kriging_variogram_exponential(h, nugget, range, sill, A) {
	return nugget + (sill - nugget) / range * (1.0 - Math.exp(-(1.0 / A) * (h / range)));
};
var kriging_variogram_spherical = function kriging_variogram_spherical(h, nugget, range, sill, A) {
	if (h > range) return nugget + (sill - nugget) / range;
	return nugget + (sill - nugget) / range * (1.5 * (h / range) - 0.5 * Math.pow(h / range, 3));
};

// Train using gaussian processes with bayesian priors
kriging.train = function (t, x, y, model, sigma2, alpha) {
	var variogram = {
		t: t,
		x: x,
		y: y,
		nugget: 0.0,
		range: 0.0,
		sill: 0.0,
		A: 1 / 3,
		n: 0
	};
	switch (model) {
		case "gaussian":
			variogram.model = kriging_variogram_gaussian;
			break;
		case "exponential":
			variogram.model = kriging_variogram_exponential;
			break;
		case "spherical":
			variogram.model = kriging_variogram_spherical;
			break;
	}

	// Lag distance/semivariance
	var i = void 0,
	    j = void 0,
	    k = void 0,
	    l = void 0,
	    n = t.length;
	var distance = Array((n * n - n) / 2);
	for (i = 0, k = 0; i < n; i++) {
		for (j = 0; j < i; j++, k++) {
			distance[k] = Array(2);
			distance[k][0] = Math.pow(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2), 0.5);
			distance[k][1] = Math.abs(t[i] - t[j]);
		}
	}distance.sort(function (a, b) {
		return a[0] - b[0];
	});
	variogram.range = distance[(n * n - n) / 2 - 1][0];

	// Bin lag distance
	var lags = (n * n - n) / 2 > 30 ? 30 : (n * n - n) / 2;
	var tolerance = variogram.range / lags;
	var lag = [0].rep(lags);
	var semi = [0].rep(lags);
	if (lags < 30) {
		for (l = 0; l < lags; l++) {
			lag[l] = distance[l][0];
			semi[l] = distance[l][1];
		}
	} else {
		for (i = 0, j = 0, k = 0, l = 0; i < lags && j < (n * n - n) / 2; i++, k = 0) {
			while (distance[j][0] <= (i + 1) * tolerance) {
				lag[l] += distance[j][0];
				semi[l] += distance[j][1];
				j++;k++;
				if (j >= (n * n - n) / 2) break;
			}
			if (k > 0) {
				lag[l] /= k;
				semi[l] /= k;
				l++;
			}
		}
		if (l < 2) return variogram; // Error: Not enough points
	}

	// Feature transformation
	n = l;
	variogram.range = lag[n - 1] - lag[0];
	var X = [1].rep(2 * n);
	var Y = Array(n);
	var A = variogram.A;
	for (i = 0; i < n; i++) {
		switch (model) {
			case "gaussian":
				X[i * 2 + 1] = 1.0 - Math.exp(-(1.0 / A) * Math.pow(lag[i] / variogram.range, 2));
				break;
			case "exponential":
				X[i * 2 + 1] = 1.0 - Math.exp(-(1.0 / A) * lag[i] / variogram.range);
				break;
			case "spherical":
				X[i * 2 + 1] = 1.5 * (lag[i] / variogram.range) - 0.5 * Math.pow(lag[i] / variogram.range, 3);
				break;
		}
		Y[i] = semi[i];
	}

	// Least squares
	var Xt = kriging_matrix_transpose(X, n, 2);
	var Z = kriging_matrix_multiply(Xt, X, 2, n, 2);
	Z = kriging_matrix_add(Z, kriging_matrix_diag(1 / alpha, 2), 2, 2);
	var cloneZ = Z.slice(0);
	if (kriging_matrix_chol(Z, 2)) kriging_matrix_chol2inv(Z, 2);else {
		kriging_matrix_solve(cloneZ, 2);
		Z = cloneZ;
	}
	var W = kriging_matrix_multiply(kriging_matrix_multiply(Z, Xt, 2, 2, n), Y, 2, n, 1);

	// Variogram parameters
	variogram.nugget = W[0];
	variogram.sill = W[1] * variogram.range + variogram.nugget;
	variogram.n = x.length;

	// Gram matrix with prior
	n = x.length;
	var K = Array(n * n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			K[i * n + j] = variogram.model(Math.pow(Math.pow(x[i] - x[j], 2) + Math.pow(y[i] - y[j], 2), 0.5), variogram.nugget, variogram.range, variogram.sill, variogram.A);
			K[j * n + i] = K[i * n + j];
		}
		K[i * n + i] = variogram.model(0, variogram.nugget, variogram.range, variogram.sill, variogram.A);
	}

	// Inverse penalized Gram matrix projected to target vector
	var C = kriging_matrix_add(K, kriging_matrix_diag(sigma2, n), n, n);
	var cloneC = C.slice(0);
	if (kriging_matrix_chol(C, n)) kriging_matrix_chol2inv(C, n);else {
		kriging_matrix_solve(cloneC, n);
		C = cloneC;
	}

	// Copy unprojected inverted matrix as K
	var K_C = C.slice(0);
	var M = kriging_matrix_multiply(C, t, n, n, 1);
	variogram.K = K_C;
	variogram.M = M;

	return variogram;
};

// Model prediction
kriging.predict = function (x, y, variogram) {
	var i = void 0,
	    k = Array(variogram.n);
	for (i = 0; i < variogram.n; i++) {
		k[i] = variogram.model(Math.pow(Math.pow(x - variogram.x[i], 2) + Math.pow(y - variogram.y[i], 2), 0.5), variogram.nugget, variogram.range, variogram.sill, variogram.A);
	}return kriging_matrix_multiply(k, variogram.M, 1, variogram.n, 1)[0];
};
kriging.variance = function (x, y, variogram) {
	var i = void 0,
	    k = Array(variogram.n);
	for (i = 0; i < variogram.n; i++) {
		k[i] = variogram.model(Math.pow(Math.pow(x - variogram.x[i], 2) + Math.pow(y - variogram.y[i], 2), 0.5), variogram.nugget, variogram.range, variogram.sill, variogram.A);
	}return variogram.model(0, variogram.nugget, variogram.range, variogram.sill, variogram.A) + kriging_matrix_multiply(kriging_matrix_multiply(k, variogram.K, 1, variogram.n, variogram.n), k, 1, variogram.n, 1)[0];
};

// Gridded matrices or contour paths
kriging.grid = function (polygons, variogram, width) {
	var i = void 0,
	    j = void 0,
	    k = void 0,
	    n = polygons.length;
	if (n == 0) return;

	// Boundaries of polygons space
	var xlim = [polygons[0][0][0], polygons[0][0][0]];
	var ylim = [polygons[0][0][1], polygons[0][0][1]];
	for (i = 0; i < n; i++) {
		// Polygons
		for (j = 0; j < polygons[i].length; j++) {
			// Vertices
			if (polygons[i][j][0] < xlim[0]) xlim[0] = polygons[i][j][0];
			if (polygons[i][j][0] > xlim[1]) xlim[1] = polygons[i][j][0];
			if (polygons[i][j][1] < ylim[0]) ylim[0] = polygons[i][j][1];
			if (polygons[i][j][1] > ylim[1]) ylim[1] = polygons[i][j][1];
		}
	} // Alloc for O(n^2) space
	var xtarget = void 0,
	    ytarget = void 0;
	var a = Array(2),
	    b = Array(2);
	var lxlim = Array(2); // Local dimensions
	var lylim = Array(2); // Local dimensions
	var x = Math.ceil((xlim[1] - xlim[0]) / width);
	var y = Math.ceil((ylim[1] - ylim[0]) / width);

	var A = Array(x + 1);
	for (i = 0; i <= x; i++) {
		A[i] = Array(y + 1);
	}for (i = 0; i < n; i++) {
		// Range for polygons[i]
		lxlim[0] = polygons[i][0][0];
		lxlim[1] = lxlim[0];
		lylim[0] = polygons[i][0][1];
		lylim[1] = lylim[0];
		for (j = 1; j < polygons[i].length; j++) {
			// Vertices
			if (polygons[i][j][0] < lxlim[0]) lxlim[0] = polygons[i][j][0];
			if (polygons[i][j][0] > lxlim[1]) lxlim[1] = polygons[i][j][0];
			if (polygons[i][j][1] < lylim[0]) lylim[0] = polygons[i][j][1];
			if (polygons[i][j][1] > lylim[1]) lylim[1] = polygons[i][j][1];
		}

		// Loop through polygon subspace
		a[0] = Math.floor((lxlim[0] - (lxlim[0] - xlim[0]) % width - xlim[0]) / width);
		a[1] = Math.ceil((lxlim[1] - (lxlim[1] - xlim[1]) % width - xlim[0]) / width);
		b[0] = Math.floor((lylim[0] - (lylim[0] - ylim[0]) % width - ylim[0]) / width);
		b[1] = Math.ceil((lylim[1] - (lylim[1] - ylim[1]) % width - ylim[0]) / width);
		for (j = a[0]; j <= a[1]; j++) {
			for (k = b[0]; k <= b[1]; k++) {
				xtarget = xlim[0] + j * width;
				ytarget = ylim[0] + k * width;
				if (polygons[i].pip(xtarget, ytarget)) A[j][k] = kriging.predict(xtarget, ytarget, variogram);
			}
		}
	}
	A.xlim = xlim;
	A.ylim = ylim;
	A.zlim = [variogram.t.min(), variogram.t.max()];
	A.width = width;
	return A;
};
kriging.contour = function (value, polygons, variogram) {
	return null;
};

// Plotting on the DOM
kriging.plot = function (canvas, grid, xlim, ylim, colors) {
	// Clear screen 
	var ctx = canvas.getContext("2d");
	ctx.clearRect(0, 0, canvas.width, canvas.height);

	// Starting boundaries
	var range = [xlim[1] - xlim[0], ylim[1] - ylim[0], grid.zlim[1] - grid.zlim[0]];
	var i = void 0,
	    j = void 0,
	    x = void 0,
	    y = void 0,
	    z = void 0;
	var n = grid.length;
	var m = grid[0].length;
	var wx = Math.ceil(grid.width * canvas.width / (xlim[1] - xlim[0]));
	var wy = Math.ceil(grid.width * canvas.height / (ylim[1] - ylim[0]));
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			if (grid[i][j] == undefined) continue;
			x = canvas.width * (i * grid.width + grid.xlim[0] - xlim[0]) / range[0];
			y = canvas.height * (1 - (j * grid.width + grid.ylim[0] - ylim[0]) / range[1]);
			z = (grid[i][j] - grid.zlim[0]) / range[2];
			if (z < 0.0) z = 0.0;
			if (z > 1.0) z = 1.0;

			ctx.fillStyle = colors[Math.floor((colors.length - 1) * z)];
			ctx.fillRect(Math.round(x - wx / 2), Math.round(y - wy / 2), wx, wy);
		}
	}
};

var addd = function addd(a, b) {
	return add(a, b);
};
var K = {
	addd: addd,
	kriging: kriging
};

function _defaults(obj, defaults) { var keys = Object.getOwnPropertyNames(defaults); for (var i = 0; i < keys.length; i++) { var key = keys[i]; var value = Object.getOwnPropertyDescriptor(defaults, key); if (value && value.configurable && obj[key] === undefined) { Object.defineProperty(obj, key, value); } } return obj; }

function _classCallCheck(instance, Constructor) { if (!(instance instanceof Constructor)) { throw new TypeError("Cannot call a class as a function"); } }

function _possibleConstructorReturn(self, call) { if (!self) { throw new ReferenceError("this hasn't been initialised - super() hasn't been called"); } return call && (typeof call === "object" || typeof call === "function") ? call : self; }

function _inherits(subClass, superClass) { if (typeof superClass !== "function" && superClass !== null) { throw new TypeError("Super expression must either be null or a function, not " + typeof superClass); } subClass.prototype = Object.create(superClass && superClass.prototype, { constructor: { value: subClass, enumerable: false, writable: true, configurable: true } }); if (superClass) Object.setPrototypeOf ? Object.setPrototypeOf(subClass, superClass) : _defaults(subClass, superClass); }

/*
width is the pixel cell's width
model is the kriging render method,mainly includ 'Gaussian','Exponential',and 'Spherical'
Gaussian: k(a,b) = w[0] + w[1] * ( 1 - exp{ -( ||a-b|| / range )2 / A } )
Exponential: k(a,b) = w[0] + w[1] * ( 1 - exp{ -( ||a-b|| / range ) / A } )
Spherical: k(a,b) = w[0] + w[1] * ( 1.5 * ( ||a-b|| / range ) - 0.5 * ( ||a-b|| / range )3 )
*/
var options = {
    width: 0.001,
    model: 'exponential',
    sigma2: 0,
    alpha: 10
};

var KrigingLayer = function (_maptalks$Layer) {
    _inherits(KrigingLayer, _maptalks$Layer);

    function KrigingLayer(id, interest, options) {
        _classCallCheck(this, KrigingLayer);

        if (!Array.isArray(interest)) {
            options = interest;
            interest = null;
        }

        var _this = _possibleConstructorReturn(this, _maptalks$Layer.call(this, id, options));

        _this._interest = interest || [];
        return _this;
    }

    KrigingLayer.prototype.getData = function getData() {
        return this._interest;
    };

    KrigingLayer.prototype.setData = function setData(interest) {
        this._interest = interest || [];
        return this.redraw();
    };

    KrigingLayer.prototype.onConfig = function onConfig(conf) {
        _maptalks$Layer.prototype.onConfig.apply(this, arguments);
        for (var p in conf) {
            if (options[p]) {
                return this.redraw();
            }
        }
        return this;
    };

    KrigingLayer.prototype.setModel = function setModel(model) {
        if (model instanceof String) {
            this.options['model'] = model;
            return this.redraw();
        }
        return this;
    };

    KrigingLayer.prototype.redraw = function redraw() {
        var renderer$$1 = this._getRenderer();
        if (renderer$$1) {
            renderer$$1.clearHeatCache();
            renderer$$1.setToRedraw();
        }
        return this;
    };

    KrigingLayer.prototype.isEmpty = function isEmpty() {
        if (!this._interest.length) {
            return true;
        }
        return false;
    };

    KrigingLayer.prototype.clear = function clear() {
        this._interest = [];
        this.redraw();
        this.fire('clear');
        return this;
    };

    /**
     * Export the KrigingLayer's JSON.
     * @return {Object} layer's JSON
     */


    KrigingLayer.prototype.toJSON = function toJSON(options) {
        if (!options) {
            options = {};
        }
        var json = {
            'type': this.getJSONType(),
            'id': this.getId(),
            'options': this.config()
        };
        var data = this.getData();
        if (options['clipExtent']) {
            var clipExtent = new maptalks.Extent(options['clipExtent']);
            var r = this._getHeatRadius();
            if (r) {
                clipExtent = clipExtent._expand(r);
            }
            var clipped = [];
            for (var i = 0, len = data.length; i < len; i++) {
                if (clipExtent.contains(new maptalks.Coordinate(data[i][0], data[i][1]))) {
                    clipped.push(data[i]);
                }
            }
            json['data'] = clipped;
        } else {
            json['data'] = data;
        }

        return json;
    };

    /**
     * Reproduce a KrigingLayer from layer's JSON.
     * @param  {Object} json - layer's JSON
     * @return {maptalks.KrigingLayer}
     * @static
     * @private
     * @function
     */


    KrigingLayer.fromJSON = function fromJSON(json) {
        if (!json || json['type'] !== 'KrigingLayer') {
            return null;
        }
        return new KrigingLayer(json['id'], json['data'], json['options']);
    };

    return KrigingLayer;
}(maptalks.Layer);

KrigingLayer.mergeOptions(options);

KrigingLayer.registerJSONType('KrigingLayer');

KrigingLayer.registerRenderer('canvas', function (_maptalks$renderer$Ca) {
    _inherits(_class, _maptalks$renderer$Ca);

    function _class() {
        _classCallCheck(this, _class);

        return _possibleConstructorReturn(this, _maptalks$renderer$Ca.apply(this, arguments));
    }

    _class.prototype.needToRedraw = function needToRedraw() {
        var map = this.layer.getMap();
        if (map.isZooming()) {
            return false;
        }
        if (map.isMoving()) {
            return false;
        }
        return _maptalks$renderer$Ca.prototype.needToRedraw.call(this);
    };

    _class.prototype.draw = function draw() {
        this.prepareCanvas();
        if (!this._isInExtent()) {
            this.completeRender();
            return;
        }
        this._plot();
        this.completeRender();
    };

    _class.prototype.drawOnInteracting = function drawOnInteracting() {
        this.draw();
    };

    _class.prototype._plot = function _plot() {
        var map = this.layer.getMap();
        var width = this.layer.options['width'];
        var colors = this.layer.options['colors'];
        var regions = this.layer.options['regions'];
        var model = this.layer.options['model'];
        var sigma2 = this.layer.options['sigma2'];
        var alpha = this.layer.options['alpha'];
        var _polygons = this._handRegions(regions);
        var extent = map.getExtent();
        var data = this.layer.getData();
        var lngs = data.map(function (d) {
            return d[0];
        });
        var lats = data.map(function (d) {
            return d[1];
        });
        var values = data.map(function (d) {
            return d[2];
        });
        var variogram = K.kriging.train(values, lngs, lats, model, sigma2, alpha);
        var grid = K.kriging.grid(_polygons, variogram, width);
        K.kriging.plot(this.canvas, grid, [extent.xmin, extent.xmax], [extent.ymin, extent.ymax], colors);
    };

    _class.prototype._isInExtent = function _isInExtent() {
        var map = this.layer.getMap();
        var mapExtent = map.getExtent();
        var regions = this.layer.options['regions'];
        var regionExtent = regions.getExtent();
        if (mapExtent.intersects(regionExtent)) return true;else return false;
    };

    _class.prototype.onZoomEnd = function onZoomEnd() {
        //delete this._heatViews;
        _maptalks$renderer$Ca.prototype.onZoomEnd.apply(this, arguments);
    };

    _class.prototype.onResize = function onResize() {
        //this._interest._width  = this.canvas.width;
        //this._interest._height = this.canvas.height;
        _maptalks$renderer$Ca.prototype.onResize.apply(this, arguments);
    };

    _class.prototype.onDragRotateEnd = function onDragRotateEnd(e) {
        _maptalks$renderer$Ca.prototype.onDragRotateEnd.call(this, e);
    };

    _class.prototype.resizeCanvas = function resizeCanvas() {
        if (!this.canvas) {
            return;
        }
    };

    _class.prototype.clearCanvas = function clearCanvas() {
        if (!this.canvas) {
            return;
        }
    };

    _class.prototype._handRegions = function _handRegions(regions) {
        var _polygons = regions.getCoordinates().map(function (coords) {
            return coords.map(function (c) {
                return [c.x, c.y];
            });
        });
        return _polygons;
    };

    return _class;
}(maptalks.renderer.CanvasRenderer));

exports.KrigingLayer = KrigingLayer;

Object.defineProperty(exports, '__esModule', { value: true });

typeof console !== 'undefined' && console.log('maptalks.kriginglayer v0.1.0, requires maptalks@^0.37.0.');

})));
