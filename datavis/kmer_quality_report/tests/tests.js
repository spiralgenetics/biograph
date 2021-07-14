/* See https://qunitjs.com/cookbook/ */

test('qUnit initialized', function() {
    var value = 'hello'
    equal(value, 'hello', 'Hello qUnit!')
})

test('dataset loaded successfully', function() {
    ok(_.every(kqr.datasets, function(d) {
        return d.ready
    }), 'all datasets are "ready"')

    _.each(kqr.datasets, function(ds) {
        var name = ds.metadata.name
        ok(name, 'dataset name: ' + name)

        var records = ds.data.length
        ok(records > 0, records + ' records loaded')
    })
})

asyncTest('chart has data', function() {
    expect(1)

    setTimeout(function() {
        ok(kqr.datapointsShown, kqr.datapointsShown + ' datapoints shown')
        start()
    }, 100)
})

test('crossfilter is working', function() {
    var count = kqr.datasets[0].byX.bottom(Infinity)
        .length
    ok(count, 'crossfilter byX.bottom has ' + count + ' records')

    var count = kqr.datasets[0].byY.top(Infinity)
        .length
    ok(count, 'crossfilter byY.top has ' + count + ' records')
})

asyncTest('chart is properly scaled', function() {
    expect(2)

    setTimeout(function() {
        equal(kqr.nvchart.xDomain()[1], kqr.xScale, 'the X scale is set')
        equal(kqr.nvchart.yDomain()[1], kqr.yScale, 'the Y scale is set')
        start()
    }, 100)
})

asyncTest('min buttons work', function() {
    expect(2)

    var prev = kqr.datapointsShown

    ok(prev, prev + ' datapoints shown')

    $('#set_min_x')
        .click()
    $('#set_min_y')
        .click()

    setTimeout(function() {
        ok(kqr.datapointsShown < prev,
            'there are fewer datapoints after clicking min (' + kqr.datapointsShown + ' < ' + prev + ')')
        start()
    }, 100)

})

asyncTest('max buttons work', function() {
    expect(2)

    var prev = kqr.datapointsShown

    ok(prev, prev + ' datapoints shown')

    $('#set_max_x')
        .click()
    $('#set_max_y')
        .click()

    setTimeout(function() {
        ok(kqr.datapointsShown > prev,
            'there are more datapoints after clicking max (' + kqr.datapointsShown + ' > ' + prev + ')')
        start()
    }, 100)

})

test('sliders work', function() {
    // Refreshing the chart sets the slider values automatically
    kqr.xScale = 100
    kqr.yScale = 20000
    set_bin_size(10)

    ok($('#set_x')
        .slider('value') == kqr.xScale, 'X slider works')
    ok($('#set_y')
        .slider('value') == kqr.yScale, 'Y slider works')
    ok($('#set_detail')
        .slider('value') == 1, 'detail slider works')

    set_bin_size(1)
})

asyncTest('auto buttons work', function() {
    expect(2)

    var prev = kqr.datapointsShown

    ok(prev, prev + ' datapoints shown')

    $('#set_ideal_x')
        .click()
    $('#set_ideal_y')
        .click()

    setTimeout(function() {
        ok(kqr.datapointsShown < prev,
            'there are fewer datapoints after clicking auto (' + kqr.datapointsShown + ' < ' + prev + ')')
        start()
    }, 100)

})