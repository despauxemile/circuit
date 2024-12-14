## circuit
exploring new ways to reinvent the wheel

### changelog
#### now
- modified nodal analysis
- dc operating point solving (R, C, L, Vs, Is)

> example circuit ([falstad.com](https://falstad.com/circuit/circuitjs.html)):
> 
> ![falstad circuit example](/res/circuit-ex1.png)
> 
> [edit live](https://falstad.com/circuit/circuitjs.html?ctz=CQAgjCAMB0l3BWcMBMcUHYMGZIA4UA2ATmIxABYIkFIQEBTAWjDACgA3cFPETFbrzAUKUMRTpI606AjYAnQeBF8EhZaLphIbAMar1w0SjV8MAzbDgQU0CqQePiJSISlWIOxWB5mBJ9X4xbQUDDRBsQiEVLR0AGwio8MjeIK0oaB9aeBzIVlIoUIC-RNTzYJ0AcyUgn1S8Xmk2FKUjJWdwMQAVAAUAHQBndhbituKOzxBewZRmpJa2lonu-oGKOd4R8pbsNE66aYG5EdMg4t30g9XsNiA)
>
> program input:
> ```rs
> let elts: Vec<Element> = Vec::from([
>   VSource(0, 1, 5.0).into(),
>   Resistor(1, 2, 10).into(),
>   Capacitor::new(2, 3, 0.00001).into(),
>   Resistor(3, 0, 10).into(),
>   Resistor(2, 4, 10).into(),
>   Inductor::new(4, 5, 1.0).into(),
>   Resistor(5, 3, 10).into(),
> ]);
> ```
> program output:
> ```bash
> $ ./circuit
> node voltages: [5.0, 3.75, 1.25, 2.5, 2.5]
> ```
