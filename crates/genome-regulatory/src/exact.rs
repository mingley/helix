//! Exact source-decimal reductions for outcome ranks, never prediction rounding.
use anyhow::{ensure, Context, Result};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, ToPrimitive};

pub fn decimal(source: &str) -> Result<BigRational> {
    let source = source.trim();
    ensure!(
        !source.is_empty() && source.len() <= 1024,
        "invalid exact source score length"
    );
    let (mantissa, exponent) = if let Some(index) = source.find(['e', 'E']) {
        (
            &source[..index],
            source[index + 1..]
                .parse::<i32>()
                .context("invalid source decimal exponent")?,
        )
    } else {
        (source, 0)
    };
    ensure!(
        exponent.unsigned_abs() <= 1000,
        "exact source exponent exceeds budget"
    );
    let (negative, mantissa) = if let Some(value) = mantissa.strip_prefix('-') {
        (true, value)
    } else {
        (false, mantissa.strip_prefix('+').unwrap_or(mantissa))
    };
    let parts: Vec<_> = mantissa.split('.').collect();
    ensure!(
        parts.len() <= 2 && parts.iter().all(|p| p.bytes().all(|b| b.is_ascii_digit())),
        "invalid source decimal score"
    );
    let fraction = parts.get(1).copied().unwrap_or("");
    let digits = format!("{}{}", parts[0], fraction);
    ensure!(!digits.is_empty(), "source decimal has no digits");
    let mut numerator =
        BigInt::parse_bytes(digits.as_bytes(), 10).context("invalid decimal numerator")?;
    if negative {
        numerator = -numerator;
    }
    let scale = fraction.len() as i32 - exponent;
    let denominator = if scale >= 0 {
        BigInt::from(10u8).pow(scale as u32)
    } else {
        numerator *= BigInt::from(10u8).pow((-scale) as u32);
        BigInt::one()
    };
    Ok(BigRational::new(numerator, denominator))
}

pub fn mean(values: &[BigRational]) -> Result<BigRational> {
    ensure!(!values.is_empty(), "empty exact source mean");
    let sum: BigRational = values.iter().cloned().sum();
    Ok(sum / BigInt::from(values.len()))
}

pub fn to_f64(value: &BigRational) -> Result<f64> {
    let number = value
        .to_f64()
        .context("exact source mean cannot be represented")?;
    ensure!(
        number.is_finite(),
        "exact source mean overflows floating metrics"
    );
    Ok(number)
}

pub fn ranks(values: &[BigRational]) -> Vec<f64> {
    let mut order: Vec<_> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].cmp(&values[b]));
    let mut ranks = vec![0.0; values.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && values[order[start]] == values[order[end]] {
            end += 1;
        }
        let rank = (start + end + 1) as f64 / 2.0;
        for &index in &order[start..end] {
            ranks[index] = rank;
        }
        start = end;
    }
    ranks
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Zero;
    #[test]
    fn source_decimals_preserve_exact_position_mean_ties() {
        let one = mean(&[
            decimal("0.01").unwrap(),
            decimal("0.02").unwrap(),
            decimal("0.03").unwrap(),
        ])
        .unwrap();
        let two = mean(&[decimal("0.02").unwrap()]).unwrap();
        assert_eq!(one, two);
        assert_eq!(
            ranks(&[one, two, decimal("0.03").unwrap()]),
            vec![1.5, 1.5, 3.0]
        );
        assert_eq!(decimal("-1.2e-2").unwrap(), decimal("-0.012").unwrap());
        assert!(decimal("nan").is_err());
        assert!(decimal("1e100000").is_err());
        assert!(decimal("-0.0").unwrap().is_zero());
    }
}
