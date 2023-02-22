from pytweezer.t_matrix import TMatrix, TMatrixMie




'''
def test_power_matrix():
    T = TMatrix.simple(...
    'sphere', 1.0, 'wavelength0', 1.0, ...
    'index_medium', 1.0, 'index_particle', 1.2);

    testCase.TestData.tol = 1.0e-6;
    T = testCase.TestData.T;
    checkWithTol(testCase, abs(1.0 - sum(abs(T.total.data).^2, 2)), ...
        'Power not conserved');
    Ttotal = testCase.TestData.T.total;
    Ttotal.Nmax = Ttotal.Nmax + 5;
    checkWithTol(testCase, abs(1.0 - sum(abs(Ttotal.data).^2, 2)), ...
        'Resizing Ttotal doesnt conserve power');
    Tscat = testCase.TestData.T.scattered;
    Tscat.Nmax = Tscat.Nmax + 5;
    checkWithTol(testCase, abs(1.0 - sum(abs(Tscat.total.data).^2, 2)), ...
        'Resizing Ttotal doesnt conserve power');


def testConvert(testCase)
    Ttotal = testCase.TestData.T.total;
    Tscat = testCase.TestData.T.scattered;

    assert(strcmpi(Ttotal.type, 'total'), 'Conversion to total');
    assert(strcmpi(Tscat.type, 'scattered'), ...
        'Conversion to scattered');

    % Check conversions back
    assert(strcmpi(Tscat.total.type, 'total'), ...
        'Conversion back to total');
    assert(strcmpi(Ttotal.scattered.type, 'scattered'), ...
        'Conversion back to scattered');

    % Check string convert
    import matlab.unittest.constraints.Matches;
    Ttotal = testCase.TestData.T;
    Ttotal.type = 'total';
    testCase.verifyThat(Ttotal.type, Matches('total'), ...
    'Incorrect type with string conversion');

def testResizing(testCase):
    T = testCase.TestData.T;
    Tnew1 = T;

    Tnew1.Nmax = Tnew1.Nmax + 5;
    assert(all(Tnew1.Nmax == T.Nmax + 5), ...
        'Faild to increase Nmax with vector size');
    assert(all(size(Tnew1.data) > size(T.data)), ...
        'Tmatrix size not actually increased (vector input)');

    Tnew2 = Tnew1;
    Tnew2.Nmax = T.Nmax;
    assert(all(Tnew2.Nmax == T.Nmax), ...
        'Failed to decrease Nmax with vector size');
    assert(all(size(Tnew2.data) == size(T.data)), ...
        'Tmatrix size not actually decreased (vector input)');

    Tnew1 = T;
    Tnew1.Nmax = T.Nmax(1) + 5;
    assert(all(Tnew1.Nmax == T.Nmax + 5), ...
        'Faild to increase Nmax (scalar input)');
    assert(all(size(Tnew1.data) > size(T.data)), ...
        'Tmatrix size not actually increased (scalar input)');

    Tnew1 = T;
    Tnew1.Nmax = [T.Nmax(1) + 5, T.Nmax(2)];
    assert(all(Tnew1.Nmax == [T.Nmax(1) + 5, T.Nmax(2)]), ...
        'Faild to increase Nmax (uneven input)');
    assert(size(Tnew1.data, 1) > size(T.data, 1) ...
        && size(Tnew1.data, 2) == size(T.data, 2), ...
        'Tmatrix size not increased correctly (uneven input)');

    Tnew1 = T;
    Tnew1.Nmax(1) = T.Nmax(1) + 5;
    assert(all(Tnew1.Nmax == [T.Nmax(1) + 5, T.Nmax(2)]), ...
        'Faild to increase Nmax (index input)');
    assert(size(Tnew1.data, 1) > size(T.data, 1) ...
        && size(Tnew1.data, 2) == size(T.data, 2), ...
        'Tmatrix size not increased correctly (index input)');

def testShrinkPowerWarning(testCase):
    T = testCase.TestData.T;
    testCase.verifyWarning(@() T.set_Nmax(1), ...
        'ott:Tmatrix:setNmax:truncation');

    T = testCase.TestData.T;
    testCase.verifyWarningFree(@() T.set_Nmax(1, 'powerloss', 'ignore'), ...
        'Ignore argument had no effect');

def testRealImagFunctions(testCase):

  T = testCase.TestData.T;

  T2 = real(T);
  T3 = imag(T);

  T4 = T.real();
  T5 = T.imag();


'''