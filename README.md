# maptalks.KrigingLayer

[![CircleCI](https://circleci.com/gh/liubgithub/maptalks.KrigingLayer/tree/master.svg?style=shield)](https://circleci.com/gh/liubgithub/maptalks.KrigingLayer)
[![NPM Version](https://img.shields.io/npm/v/maptalks.KrigingLayer.svg)](https://github.com/liubgithub/maptalks.KrigingLayer)

A maptalks Layer to draw spatial data in a way of providing spatial prediction and mapping capabilities via the ordinary kriging algorithm.

![screenshot](https://liubgithub.github.io/archives/kriging/screenshot.png)

## Examples

* [Lintong ,Xi'an's kriging map](https://liubgithub.github.io/archives/kriging/). (inspired by [oeo4b](https://github.com/oeo4b/kriging.js)).In this example,the red points are the predict sample points.
* [the world of kriging map](http://oeo4b.github.io).

## Install
  
* Install with npm: ```npm install maptalks.KrigingLayer```. 
* Download from [dist directory](https://github.com/liubgithub/maptalks.KrigingLayer/tree/master/dist).
* Use unpkg CDN: ```https://unpkg.com/maptalks.KrigingLayer/dist/maptalks.KrigingLayer.min.js```

## Usage

As a plugin, ```maptalks.KrigingLayer``` must be loaded after ```maptalks.js``` in browsers.
```html
<link rel="stylesheet" href="https://unpkg.com/maptalks/dist/maptalks.css">
<script type="text/javascript" src="https://unpkg.com/maptalks/dist/maptalks.min.js"></script>
<script type="text/javascript" src="https://unpkg.com/maptalks.KrigingLayer/dist/maptalks.KrigingLayer.min.js"></script>
<script>
    var samples = [[109.199, 34.443, 4.278], [109.224, 34.464, 92.767], [109.271, 34.464, 45.834], [109.322, 34.447, 5.220]];
    var regions =[[109.220, 34.589], [109.161, 34.648], [109.131, 34.631], [109.156, 34.526]
    , [109.170, 34.437]];
    var colors =["#00A600", "#01A600", "#03A700", "#04A700", "#05A800", "#07A800"];
    var polygon = new maptalks.Polygon(regions);
    var krigingLayer = new maptalks.KrigingLayer('kriging', samples, {
        colors: colors,
        regions: polygon
    });
</script>
```
## Supported Browsers

IE 9-11, Chrome, Firefox, other modern and mobile browsers.

## API Reference

```KrigingLayer``` is a subclass of [maptalks.Layer](http://maptalks.org/api.maptalks.org/master/Layer.html) and inherits all the methods of its parent.

### `Constructor`

```javascript
// samples's format
//[[x,y,value], [x,y,value] ..]
// symbol only supports lineWidth and lineColor
new maptalks.KrigingLayer(id, samples, options)
```

* id **String** layer id
* samples **Object[]** providing sample data for trainning method, `[[x,y,value], [x,y,value] ..]`
* options **Object** options
    * alpha **Number** is a number that correspond to the variance parameters of the gaussian process and the prior of the variogram model, respectively.A diffuse α prior is typically used (10 by default)
    * sigma2 **Number** represent σ2,it's a variance parameters of the gaussian process  (0 by default)
    * width **Number** the pixel cell's width (0.001 by default)
    * model **String** the kriging render method,mainly includ 'Gaussian','Exponential',and 'Spherical' (exponential by default)
    * colors **Object[]** the color you will render the layer
    * regions **maptlaks.Polygon** render polygons

### `getData()`

get layer's data

**Returns** `Object[]`

### `setData(data)`

set new data to the layer

* data **Object[]** new predict points

**Returns** `this`

### `setModel(model)`

set new model to the kriging render method

* model **String** render method

**Returns** `this`


## Contributing

We welcome any kind of contributions including issue reportings, pull requests, documentation corrections, feature requests and any other helps.

## Develop

The only source file is ```index.js```.

It is written in ES6, transpiled by [babel](https://babeljs.io/) and tested with [mocha](https://mochajs.org) and [expect.js](https://github.com/Automattic/expect.js).

### Scripts

* Install dependencies
```shell
$ npm install
```

* Watch source changes and generate runnable bundle repeatedly
```shell
$ npm run dev
```

* Tests
```shell
$ npm test
```

* Watch source changes and run tests repeatedly
```shell
$ npm run tdd
```

* Package and generate minified bundles to dist directory
```shell
$ npm run build
```

* Lint
```shell
$ npm run lint
```
