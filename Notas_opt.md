
# Notas de algoritmo de Otimização
# 1. Introdução
<h3> 

$\hspace{1cm}$ De forma geral, o problema tratado terá a seguinte forma: 

<br>
<center>

minimizar  $f(x) $

sujeito a $x$ $ \epsilon$ $ X$

</center>
</h3>

## 1.1 - Regressão Rigde

<h3> 

$\hspace{1cm}$ Dado os pontos $(x_1,y_1),(x_2,y_2),...,(x_N,y_N)$ queremos a melhor relação entre y e x que atenda $y  \approx x^T \beta$, em que $\beta$ são os pesos do modelo, que precisam ser encontrados.

Assim, teremos:

<center>

minimizar  $\frac{1}{N}\sum (y_i - x_i\beta) + \lambda\sum_j^d \beta_j^2 $, &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;        (1)

sujeito a $\beta$ $ \epsilon$ $ \R$


</center>

o primeiro termo acima o erro médio, e o segundo é a penalização do modelo.


$\hspace{1cm}$ Há duas maneiras de encontrar os valores dos pesos, a primeira é usando o algoritmo do gradiente descendente e a segunda é o método de Newton,


$\hspace{0.5cm}$ Gradiente descendente $\hspace{2.5cm}$ Newton

$\hspace{0.5cm}$  $\beta_{i + 1} = \beta_i - \eta \nabla_\beta f(\beta_i)$ $\hspace{0.5cm}$ $\beta_{i + 1} = \beta_i - \eta(\nabla_\beta^2 f(\beta_i))^{-1} \nabla_\beta f(\beta_i)$

Um ponto importante, é que G.D. segue uma trajetória perpendicular as linhas de nível, enquanto Newton vai direto ao mínimo.

![GD](/images/GD.png)

Contudo,newton se torna mais custoso a medida que computar $(\nabla_\beta^2 f(\beta_i))^{-1}$ para um conjunto de daods muito grande se torna difícil.

Não apenas isto, mas haveão funções que não possuiram segunda derivadas, e assim impossibilitando o metodo de newton. O exemplo disto está na regressão de Lasso.
</h3>

## 1.2 - Regressão Lasso

<h3>

 $\hspace{1cm}$ Dado os pontos $(x_1,y_1),(x_2,y_2),...,(x_N,y_N)$ queremos a melhor relação entre y e x que atenda $y  \approx x^T \beta$, em que $\beta$ são os pesos do modelo, que precisam ser encontrados.

Assim, teremos:

<center>

minimizar  $\frac{1}{N}\sum (y_i - x_i\beta) + \lambda\sum_j^d |\beta_j| $, &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;        (2)

sujeito a $\beta$ $ \epsilon$ $ \R$

</center>

$\hspace{1cm}$ o termo de penalidade não possui derivada em zero. Veremos mais a frente como lidar com isto.
</h3>

# 2. Básico de Algebra linear e cálculo

## 2.1 - Algebra Linear
<h3>

$\hspace{1cm}$ Supondo um vetor $v$, quando dizemos que tal pertence a $\R^d$, siginifica que $ v=(v_1,v_2,...,v_d)$

$\hspace{1cm}$ Supondo dois vetores no espeço bidimensional, u e v, dizemos que {u,v} spam $\R^2$, se é possível escrever um terceiro vetor w qualquer como, $w = \alpha u + \beta v$. É possivel 
generalizar isto para o espaço $\R^k$, em que agora, $w = \sum_i^k \alpha _i v_i$.

 $\hspace{1cm}$ Com isto é possível definir a independencia linear, ou seja, supondo dois vetores u e v, se não existir $\alpha,\beta$ que consiga representar um em função do  outro, dizemos que eles são linearmente independentes. Uma forma matemática de descrever isto é

$\alpha u + \beta v = 0 \rightarrow  u = - (\beta/\alpha)v$

portanto, não pode haver qualquer alpha igual a zero!!!

Assim, se $(v_1,v_2,...,v_d)$ são linearmente independentes e spam $\R^d$, então eles formam uma base.

</h3>

## 2.1.1 - Operadores

<h3>

*  produto interno

$$
u \cdot v = <u,v> = \sum_i u_iv_i;
$$

isto gera um scalar.

<br>

* Produto externo

$$
u @ v^T = (u_1 u_2 ... u_N) @ (v_1 v_2 ... v_N)^T;
$$

isto gera uma matrix NxN.

<br>

* Cauchy Schwartz

$$

<u,v> \leq ||u||_2 \cdot ||v||_2 ;

$$

<br>

* Angulo

$$

<u,v> =\frac{cos(\theta)}{||u||_2 \cdot ||v||_2 } ;

$$

<br>

* Projeção: supondo dois vetores u e v, a projeção de u usando v

$$

u = \frac{<u,v>}{||u||_2} \cdot u;

$$

<br>

* Projeção em sub espaço

$$
Proj(v) = <v,u_1> u_1 + <v,u_2> u_2 + ...+ <v,u_d> u_d;
$$

<br>
* L-p norma

$$
||x||_p = (\sum(|x_i|^p))^{1/p};
$$

<br>

* Matriz rank, determinante, espaço vazio e alcance (range): Supondo uma amtriz M(m,n)

  * Null space (M): ${[x \epsilon \R^n : Mx = 0]}$;
  * Range (M): ${[y \epsilon \R^m : Mx = y, \text{for x in }\R^n]}$;
  * Rank (M): (n,m);
  * Se det(M) = 0, então as colunas de M não são independentes ou as linhas não são.
<br>

* Subespaço: V é um subespaço se ele le fechado sob adição e scaling.

<br>

* Complemento ortogonal do Subespaço: para um subespaço V $$ V^O = [w: <w,v> = 0,\forall\ v\ \epsilon \ V] $$

<br>

* Affine subespaço: Subsespaço descolado por umvetor fixo

<br>

* Matrizes simetricas: Supondo uma  Matriz M simetrica

    * M é (n,n);
    * $M = M^T$;
    * todos autovaroes são perpendiculares e todos os autovalores são reais;
    * Podemos escrever: $M = \sum\lambda_iv_iv_i^T$, sendo v_i os autovetores e lambda_i os autovalores;
<br>

* Matrizes semidefinidas e positivas: Se M é simetrica e todos os autovalores são positivos, então M é dita semi definida postiva


</h3>


## 2.2 - Cálculo
<h3>

* gradiente:

$$

\nabla f(x) = \begin{pmatrix}
\frac{\partial f(x)}{\partial x_1} \\
\frac{\partial f(x)}{\partial x_2}  \\
... \\
\frac{\partial f(x)}{\partial x_n} 
\end{pmatrix}

$$

<br>

* Hessian:

$$

\nabla^2 f(x) = \begin{pmatrix}
\frac{\partial^2 f(x)}{\partial x_1x_1} & ... & ... &  \frac{\partial^2 f(x)}{\partial x_1x_n}\\
\frac{\partial^2 f(x)}{\partial x_2x_1} & ... & ... & \frac{\partial^2 f(x)}{\partial x_2x_n}  \\
...  & ... & ...  & ... \\
\frac{\partial^2 f(x)}{\partial x_nx_1} & ... & ... & \frac{\partial^2 f(x)}{\partial x_nx_n}
\end{pmatrix}

$$

</h3>

# 3. Conjuntos convexos

## 3.1 - Definição
<h3>

$\hspace{1cm}$ É dificil dar uma definição formal sobre conjuntos conveos, porém, visualmente, conjuntos convexos são aqueles nos quais é possível traçar uma reta de um ponto a outro dentro do conjunto, e todos os pontos contidos nesta reta pertecem ao conjunto. 

$\hspace{1cm}$ Isto está ilustrado na figura abaixo.

![GD](/images/convex.png)

$\hspace{1cm}$ matematicamente, seja C um conjunto convexo, se $\forall \ x,y \ \epsilon \ C$, $\forall \ \lambda \ \epsilon \ [0,1]$, então o conjunto C será convex se, $(\lambda x + (1-\lambda)y) \ \epsilon \  C$








