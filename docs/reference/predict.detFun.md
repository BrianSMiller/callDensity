# predict method for detFun objects

Ensures that base R's stats::predict() works correctly on objects that
have been tagged with the "detFun" class by fitDetFun(). Without this,
S4 dispatch for vglm objects fails because R sees "detFun" first and
finds no registered S4 method for that class.

## Usage

``` r
# S3 method for class 'detFun'
predict(object, ...)
```

## Arguments

- object:

  A fitted detFun object from fitDetFun().

- ...:

  Passed to the underlying model's predict method.

## Details

This method strips "detFun" from the class vector and redispatches. For
vglm (S4), it calls VGAM::predict directly to bypass S3/S4 dispatch
confusion.
